#!/usr/bin/env python

import os
import sys
import time
import json
import argparse
import subprocess as sp
import pickle
from tempfile import TemporaryDirectory

import warnings
warnings.filterwarnings("ignore")

import numpy as np
from sklearn.decomposition import PCA

import mdtraj
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

WORK_HOME = os.getenv("PREFMD_HOME")
assert WORK_HOME is not None
sys.path.insert(0, '%s/bin'%WORK_HOME)

import path
from libcommon import *

from libquat import Quaternion
from libcustom import *

def solvate_pdb(output_prefix, pdb, options):
    orient_fn = path.Path('%s.orient.pdb'%output_prefix)
    if not orient_fn.status():
        xyz = pdb.xyz[0]
        xyz -= xyz.mean(axis=0)
        #
        for i in range(2):
            pca = PCA(n_components=(i+1))
            pca.fit(xyz)
            #
            axis_0 = pca.components_[i]
            axis_1 = np.zeros(3, dtype=float)
            axis_1[i] = 1.
            #
            axis_r = np.cross(axis_0, axis_1)
            angle_r = np.arccos(np.dot(axis_0, axis_1))
            #
            q = Quaternion.from_axis_and_angle(axis_r, angle_r)
            #
            xyz = np.dot(xyz, q.rotate().T)
        #
        system_size = np.max(xyz, axis=0) - np.min(xyz, axis=0) + 2.0*(options['md']['solvate'] * 0.1)
        translate = system_size / 2.0
        xyz += translate
        pdb.xyz[0] = xyz
        #
        pdb.unitcell_vectors = (system_size * np.eye(3))[None,:]
        pdb.save(orient_fn.short())
    #
    solv_fn = path.Path('%s.solvate.pdb'%output_prefix)
    if not solv_fn.status():
        cmd = []
        cmd.append("%s/solvate.py"%EXEC_HOME)
        cmd.append(orient_fn.short())
        cmd.append(solv_fn.short())
        cmd.append("%8.5f"%options['md']['ion_conc'])
        system(cmd)
        #
        cmd = []
        cmd.append("%s/update_water_name.py"%EXEC_HOME)
        cmd.append(solv_fn.short())
        system(cmd)
        #
        #if options['ff']['use_modified_CMAP']:
        #    cmd = ['%s/resName_modified_CMAP.py'%EXEC_HOME, solv_fn.short()]
        #    output = system(cmd, verbose=verbose, stdout=True)
        #    with solv_fn.open('wt') as fout:
        #        if 'ssbond' in options:
        #            for line in options['ssbond']:
        #                fout.write("%s\n"%line)
        #        fout.write(output)
    #
    return orient_fn, solv_fn

def generate_PSF(output_prefix, solv_fn, options):
    psf_fn = path.Path("%s.psf"%output_prefix)
    crd_fn = path.Path("%s.crd"%output_prefix)
    if psf_fn.status() and crd_fn.status():
        return psf_fn, crd_fn
    #
    cmd = []
    cmd.append("%s/genPSF.py"%EXEC_HOME)
    cmd.append(solv_fn.short())
    cmd.extend(['-psf', psf_fn.short()])
    cmd.extend(['-crd', crd_fn.short()])
    cmd.append("--toppar")
    cmd.extend(options['ff']['toppar'])
    system(cmd)
    #
    return psf_fn, crd_fn

def construct_restraint(psf, pdb, force_const):
    rsr = CustomExternalForce("k0*dsq ; dsq=((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    rsr.addPerParticleParameter("x0")
    rsr.addPerParticleParameter("y0")
    rsr.addPerParticleParameter("z0")
    rsr.addPerParticleParameter('k0')
    #
    calphaIndex = []
    for i,atom in enumerate(psf.topology.atoms()):
        if atom.name == 'CA':
            calphaIndex.append(i)
    #
    k = -1
    for i,atom in enumerate(pdb.top.atoms):
        if atom.name != 'CA': continue
        #
        k += 1
        mass = atom.element.mass
        param = pdb.xyz[0,i].tolist()
        param.append(force_const * mass * kilocalories_per_mole/angstroms**2)
        rsr.addParticle(calphaIndex[k], param)
    return rsr

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
    avrg, cntr = get_input_structures(arg, options)
    #
    orient_fn, solv_fn = solvate_pdb(arg.output_prefix, cntr, options)
    #
    psf_fn, crd_fn = generate_PSF(arg.output_prefix, solv_fn, options)
    #
    final = run_md(arg.output_prefix, solv_fn, avrg, psf_fn, crd_fn, options)
    return final

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
    output.save("%s.pdb"%arg.output_prefix)

if __name__=='__main__':
    main()
