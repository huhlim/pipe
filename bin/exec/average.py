#!/usr/bin/env python

import os
import sys
import json
import argparse
from tempfile import TemporaryDirectory

import warnings
warnings.filterwarnings("ignore")

import numpy as np
from importlib import import_module

import mdtraj
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

WORK_HOME = os.getenv("PREFMD_HOME")
assert WORK_HOME is not None
sys.path.insert(0, '%s/bin'%WORK_HOME)

import path
from libcommon import *

from libcustom import *
from libmd import solvate_pdb, generate_PSF, construct_restraint

def run_md(output_prefix, solv_fn, avrg, psf_fn, crd_fn, options):
    psf = CharmmPsfFile(psf_fn.short())
    crd = CharmmCrdFile(crd_fn.short())
    #
    pdb = mdtraj.load(solv_fn.short())
    avrg = avrg.superpose(pdb, atom_indices=avrg.top.select("name CA"),\
                               ref_atom_indices=pdb.top.select("name CA"))
    #
    box = pdb.unitcell_lengths[0]
    psf.setBox(*box)
    #
    ff = CharmmParameterSet(*options['ff']['toppar'])
    platform = Platform.getPlatformByName(options['openmm']['platform'])
    #
    sys = psf.createSystem(ff, \
                           nonbondedMethod=PME, \
                           switchDistance=0.8*nanometers, \
                           nonbondedCutoff=1.0*nanometers, \
                           constraints=HBonds)
    #
    sys.addForce(construct_restraint(psf, avrg, 1.0))
    #
    if 'custom' in options['ff'] and options['ff']['custom'] is not None:
        custom_restrains = read_custom_restraint(options['ff']['custom'])
        custom_s = construct_custom_restraint(pdb, custom_restraints[1])
        for custom in custom_s:
            sys.addForce(custom)
    #
    temp = options['md']['heat'][0]
    steps = options['md']['heat'][1]
    #
    for i in range(len(options['md']['heat'][0])):
        temp = options['md']['heat'][0][i]
        steps = options['md']['heat'][1][i]
        #
        integrator = LangevinIntegrator(temp*kelvin, \
                                        options['md']['langfbeta']/picosecond, \
                                        options['md']['dyntstep']*picosecond)
        #
        simulation = Simulation(psf.topology, sys, integrator, platform)
        simulation.context.setPositions(crd.positions)
        if i == 0:
            simulation.minimizeEnergy(maxIterations=500)
            simulation.context.setVelocitiesToTemperature(temp*kelvin)
        else:
            with open(chk_fn, 'rb') as fp:
                simulation.context.loadCheckpoint(fp.read())
        simulation.step(steps)
        #
        pdb_fn = 'heat.%03d.pdb'%temp
        simulation.reporters.append(PDBReporter(pdb_fn, steps))
        simulation.step(steps)
        #
        chk_fn = '%s.heat.restart'%(output_prefix)
        with open(chk_fn, 'wb') as fout:
            fout.write(simulation.context.createCheckpoint())
        simulation = None   # have to free CUDA-related variables
    #
    final = mdtraj.load(pdb_fn)
    final = final.atom_slice(final.top.select("protein"))
    #
    output = avrg
    output.xyz = final.xyz
    #
    return output

def read_score(fn, score_fn='RWplus'):
    index = ['RWplus', 'dDFIRE', 'DFIRE'].index(score_fn)
    #
    score = []
    with fn.open() as fp:
        for line in fp:
            if line.startswith("#"): continue
            score.append(line.strip().split()[index])
    return np.array(score, dtype=float)

def read_qual(fn):
    rmsd = []
    with fn.open() as fp:
        for line in fp:
            if line.startswith("#"): continue
            rmsd.append(line.strip().split()[0])    # INDEX needed
    return np.array(rmsd, dtype=float)

def get_input_structures(arg, options):
    top = mdtraj.load(arg.top_pdb.short())
    #
    traj_s = []
    score_s = []
    rmsd_s = []
    for i,traj_fn in enumerate(arg.input_dcd_s):
        traj = mdtraj.load(traj_fn.short(), top=top)
        traj_s.append(traj)
        n_frame = len(traj)
        #
        if options['rule'][0] in ['score', 'casp12']:
            score = read_score(arg.input_score_s[i], score_fn=options['rule'][1][0])
            assert n_frame == score.shape[0]
            score_s.append(score)
        #
        if options['rule'][0] in ['casp12']:
            rmsd = read_qual(arg.input_qual_s[i])
            assert n_frame == rmsd.shape[0]
            rmsd_s.append(rmsd)
    #
    traj_s = mdtraj.join(traj_s, check_topology=False)
    traj_s.superpose(top, atom_indices=top.topology.select("name CA"))
    #
    if options['rule'][0] == 'score':
        score_s = np.concatenate(score_s)
        score_cutoff = np.percentile(score_s, options['rule'][1][1])
        frame = (score_s < score_cutoff)
    #
    avrg = copy.deepcopy(top)
    avrg.xyz = np.mean(traj_s[frame].xyz, 0)
    #
    rmsd = mdtraj.rmsd(traj_s[frame], avrg)
    i_min = np.argmin(rmsd)
    cntr = traj_s[frame][i_min]
    #
    return avrg, cntr

def run(arg, options):
    if options['rule'][0] == 'cluster':
        cluster_s = import_module("libcluster").get_clusters(arg.top_pdb, arg.input_dcd_s,
                rmsd_cutoff=options['rule'][1][0], subsample=options['rule'][1][1])
        final_s = []
        for i in range(min(options['rule'][1][2], len(cluster_s))):
            output_prefix = '%s.%04d'%(arg.output_prefix, i)
            cntr = cluster_s[i][1]
            avrg = cluster_s[i][2]
            #
            orient_fn, solv_fn = solvate_pdb(output_prefix, cntr, options, False)
            psf_fn, crd_fn = generate_PSF(output_prefix, solv_fn, options, False)
            final = run_md(output_prefix, solv_fn, avrg, psf_fn, crd_fn, options)
            final_s.append((output_prefix, final))
    else:
        avrg, cntr = get_input_structures(arg, options)
        #
        orient_fn, solv_fn = solvate_pdb(arg.output_prefix, cntr, options, False)
        psf_fn, crd_fn = generate_PSF(arg.output_prefix, solv_fn, options, False)
        final = run_md(arg.output_prefix, solv_fn, avrg, psf_fn, crd_fn, options)
        final_s = [(arg.output_prefix, final)]
    return final_s

def main():
    arg = argparse.ArgumentParser(prog='average')
    arg.add_argument(dest='output_prefix')
    arg.add_argument(dest='top_pdb')
    arg.add_argument('--input', dest='input_json', required=True)
    arg.add_argument('--dcd', dest='input_dcd_s', nargs='*', default=None, required=True)
    arg.add_argument('--score', dest='input_score_s', nargs='*', default=None)
    arg.add_argument('--qual', dest='input_qual_s', nargs='*', default=None)
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    arg.top_pdb = path.Path(arg.top_pdb)
    if arg.input_dcd_s is not None:
        arg.input_dcd_s = [path.Path(fn) for fn in arg.input_dcd_s]
    if arg.input_score_s is not None:
        arg.input_score_s = [path.Path(fn) for fn in arg.input_score_s]
    if arg.input_qual_s is not None:
        arg.input_qual_s = [path.Path(fn) for fn in arg.input_qual_s]
    #
    with open(arg.input_json) as fp:
        options = json.load(fp)
    #
    cwd = os.getcwd()
    tmpdir = TemporaryDirectory(prefix='avrg.')
    os.chdir(tmpdir.name)
    #
    output = run(arg, options)
    #
    os.chdir(cwd)
    for prefix, final in output:
        final.save("%s.pdb"%prefix)

if __name__=='__main__':
    main()
