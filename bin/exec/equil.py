#!/usr/bin/env python

import os
import sys
import time
import json
import argparse
import subprocess as sp
import pickle

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

def solvate_pdb(output_prefix, input_pdb, options, verbose):
    orient_fn = path.Path('%s.orient.pdb'%output_prefix)
    if not orient_fn.status():
        pdb = mdtraj.load(input_pdb.short())
        #
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
        system(cmd, verbose=verbose)
        #
        cmd = []
        cmd.append("%s/update_water_name.py"%EXEC_HOME)
        cmd.append(solv_fn.short())
        system(cmd, verbose=verbose)
        #
        if options['ff']['use_modified_CMAP']:
            cmd = ['%s/resName_modified_CMAP.py'%EXEC_HOME, solv_fn.short()]
            output = system(cmd, verbose=verbose, stdout=True)
            with solv_fn.open('wt') as fout:
                if 'ssbond' in options:
                    for line in options['ssbond']:
                        fout.write("%s\n"%line)
                fout.write(output)
    #
    return orient_fn, solv_fn

def generate_PSF(output_prefix, solv_fn, options, verbose):
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
    system(cmd, verbose=verbose)
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

def equil_md(output_prefix, solv_fn, psf_fn, crd_fn, options, verbose):
    psf = CharmmPsfFile(psf_fn.short())
    crd = CharmmCrdFile(crd_fn.short())
    #
    pdb = mdtraj.load(solv_fn.short())
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
    sys.addForce(construct_restraint(psf, pdb, 0.5))
    #
    if 'custom' in options['ff'] and options['ff']['custom'] is not None:
        custom_restrains = read_custom_restraint(options['ff']['custom'])
        custom_s = construct_custom_restraint(pdb, custom_restraints[1])
        for custom in custom_s:
            sys.addForce(custom)
    #
    steps_left = options['md']['equil'][0]
    steps_heat = options['md']['heat'][0]
    temp = options['md']['heat'][1]
    temp_incr = options['md']['heat'][2]
    #
    i = 0
    while (temp < options['md']['dyntemp']):
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
        simulation.reporters.append(StateDataReporter('%s.heat.%03d.log'%(output_prefix, temp), 500, step=True, \
            time=True, kineticEnergy=True, potentialEnergy=True, temperature=True, progress=True, \
            remainingTime=True, speed=True, volume=True, totalSteps=steps_heat, separator='\t'))
        #
        simulation.step(steps_heat)
        #
        chk_fn = '%s.heat.restart'%(output_prefix)
        with open(chk_fn, 'wb') as fout:
            fout.write(simulation.context.createCheckpoint())
        simulation = None   # have to free CUDA-related variables
        temp += temp_incr ; i += 1
        steps_left -= steps_heat
    #
    sys.addForce(MonteCarloBarostat(1.0*bar, options['md']['dyntemp']*kelvin))
    integrator = LangevinIntegrator(options['md']['dyntemp']*kelvin, \
                                    options['md']['langfbeta']/picosecond, \
                                    options['md']['dyntstep']*picosecond)
    #
    simulation = Simulation(psf.topology, sys, integrator, platform)
    simulation.context.setPositions(crd.positions)
    with open(chk_fn, 'rb') as fp:
        simulation.context.loadCheckpoint(fp.read())

    simulation.reporters.append(StateDataReporter('%s.equil.log'%output_prefix, 2500, step=True, \
        time=True, kineticEnergy=True, potentialEnergy=True, temperature=True, progress=True, \
        remainingTime=True, speed=True, volume=True, totalSteps=steps_left, separator='\t'))
        #
    simulation.step(steps_left)
    #
    chk_fn = '%s.equil.restart'%(output_prefix)
    #with open(chk_fn, 'wb') as fout:
    #    fout.write(simulation.context.createCheckpoint())

    state = simulation.context.getState(getPositions=True,\
                                        getVelocities=True,\
                                        getForces=True,\
                                        getEnergy=True, \
                                        enforcePeriodicBox=True)

    with open("%s.pkl"%chk_fn, 'wb') as fout:
        pickle.dump(state, fout)
    #
    boxinfo = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(nanometer)
    #with open("boxsize", 'wt') as fout:
    #    fout.write("%13.8f %13.8f %13.8f"%(boxinfo[0][0], boxinfo[1][1], boxinfo[2][2]))
    #
    equil_pdb_fn = path.Path('%s.equil.pdb'%output_prefix)
    pdb.xyz = state.getPositions(asNumpy=True).value_in_unit(nanometer)[None,:]
    #pdb.unitcell_vectors = boxinfo[None,:]/10.0
    pdb.unitcell_vectors = boxinfo[None,:]
    pdb.save(equil_pdb_fn.short())
    #
    cmd = []
    cmd.append("%s/update_water_name.py"%EXEC_HOME)
    cmd.append(equil_pdb_fn.short())
    system(cmd, verbose=verbose)
    #
    if options['ff']['use_modified_CMAP']:
        cmd = ['%s/resName_modified_CMAP.py'%EXEC_HOME, equil_pdb_fn.short()]
        output = system(cmd, verbose=verbose, stdout=True)
        with equil_pdb_fn.open('wt') as fout:
            if 'ssbond' in options:
                for line in options['ssbond']:
                    fout.write("%s\n"%line)
            fout.write(output)

def run(input_pdb, output_prefix, options, verbose, nonstd):
    tempfile_s = []
    #
    orient_fn, solv_fn = solvate_pdb(output_prefix, input_pdb, options, verbose)
    tempfile_s.extend([orient_fn, solv_fn])
    #
    psf_fn, crd_fn = generate_PSF(output_prefix, solv_fn, options, verbose)
    tempfile_s.append(crd_fn)
    #
    equil_md(output_prefix, solv_fn, psf_fn, crd_fn, options, verbose)

def main():
    arg = argparse.ArgumentParser(prog='equil')
    arg.add_argument(dest='output_prefix')
    arg.add_argument(dest='input_pdb')
    arg.add_argument('--input', dest='input_json', required=True)
    arg.add_argument('--toppar', dest='toppar', nargs='*', default=None)
    arg.add_argument('--custom', dest='custom_file', default=None)
    arg.add_argument('--temp', dest='temp', default=None, type=float)
    arg.add_argument('--non_standard', dest='non_standard', default=False, action='store_true')
    arg.add_argument('-v', '--verbose', default=False, action='store_true')
    arg.add_argument('--keep', dest='keep', action='store_true', default=False,\
            help='set temporary file mode (default=False)')
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    #
    input_pdb = path.Path(arg.input_pdb)
    #
    with open(arg.input_json) as fp:
        options = json.load(fp)
    if arg.toppar is not None:
        options['ff']['toppar']= arg.toppar
    if arg.custom_file is not None:
        options['ff']['custom'] = arg.custom_file
    if arg.temp is not None:
        options['md']['dyntemp'] = arg.temp
    #
    run(input_pdb, arg.output_prefix, options, arg.verbose, arg.non_standard)

if __name__=='__main__':
    main()
