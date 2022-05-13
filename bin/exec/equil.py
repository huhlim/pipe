#!/usr/bin/env python

import os
import sys
import json
import argparse
import pickle
import numpy as np

import warnings
warnings.filterwarnings("ignore")

import mdtraj
try:
    from openmm.unit import *
    from openmm.openmm import *
    from openmm.app import *
except:
    from simtk.unit import *
    from simtk.openmm import *
    from simtk.openmm.app import *

WORK_HOME = os.getenv("PREFMD_HOME")
assert WORK_HOME is not None
sys.path.insert(0, '%s/bin'%WORK_HOME)

import path
from libcommon import *
from libligand import read_ligand_json, add_ligand, update_ligand_name, get_ligand_restratint

from libcustom import *
from libmd import solvate_pdb, generate_PSF, construct_restraint, construct_water_restraint, \
        construct_ligand_restraint, update_residue_name

def equil_md(output_prefix, solv_fn, psf_fn, crd_fn, options, verbose):
    psf = CharmmPsfFile(psf_fn.short())
    crd = CharmmCrdFile(crd_fn.short())
    #
    pdb = mdtraj.load(solv_fn.short())
    #
    if options['md']['solvate'] > 0.0 and options['md'].get("orient", True):  # solvated via libmd.solvate_pdb
        box = np.array(crd.positions.value_in_unit(nanometers), dtype=float)
        boxsize = np.max(box, 0) - np.min(box, 0)
    else:   # solvated prior to the script
        boxsize = pdb.unitcell_lengths[0]
    psf.setBox(*boxsize)
    #
    ff_file_s = options['ff']['toppar']
    if 'ligand' in options and len(options['ligand']['str_fn_s']) > 0:
        ff_file_s.extend(options['ff']['cgenff'])
        ff_file_s.extend([fn.short() for fn in options['ligand']['str_fn_s']])

    ff = CharmmParameterSet(*ff_file_s)
    properties = {}
    cuda_devices = os.getenv("CUDA_VISIBLE_DEVICES")
    if cuda_devices is not None:
        n_gpu = len(cuda_devices.split(","))
    else:
        n_gpu = 0
        options['openmm']['platform'] = 'CPU'
    if n_gpu > 1:
        properties['DeviceIndex'] = ",".join(['%d'%index for index in range(n_gpu)])
    platform = Platform.getPlatformByName(options['openmm']['platform'])
    if options['openmm']['platform'] == 'CPU':
        properties = {}
    #
    sys = psf.createSystem(ff, \
                           nonbondedMethod=PME, \
                           switchDistance=0.8*nanometers, \
                           nonbondedCutoff=1.0*nanometers, \
                           constraints=HBonds)
    #
    sys.addForce(construct_restraint(psf, pdb, 0.5))
    if options['md'].get("solvate_cryst", None) is not None:
        n_cryst_water = options['md'].get("n_cryst_water", -1)
        if n_cryst_water < 0:
            raise NotImplementedError
        sys.addForce(construct_water_restraint(psf, pdb, n_cryst_water, 0.5))
    #
    if 'custom' in options['ff'] and options['ff']['custom'] is not None:
        custom_restraints = read_custom_restraint(options['ff']['custom'])
        custom_s = construct_custom_restraint(pdb, custom_restraints[1])
        for custom in custom_s:
            sys.addForce(custom)
    #
    if 'ligand' in options['ff']:
        ligand_restraints = construct_ligand_restraint(options['ff']['ligand'])
        sys.addForce(ligand_restraints)
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
        simulation = Simulation(psf.topology, sys, integrator, platform, properties)
        simulation.context.setPositions(crd.positions)
        if i == 0:
            #state = simulation.context.getState(getEnergy=True)
            #print (state.getPotentialEnergy())
            simulation.minimizeEnergy(maxIterations=500)
            #state = simulation.context.getState(getEnergy=True)
            #print (state.getPotentialEnergy())
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
    if np.std(boxsize) < 0.1:
        sys.addForce(MonteCarloBarostat(1.0*bar, options['md']['dyntemp']*kelvin))
    else:
        sys.addForce(MonteCarloAnisotropicBarostat([1.0*bar, 1.0*bar, 1.0*bar], \
                options['md']['dyntemp']*kelvin))
    integrator = LangevinIntegrator(options['md']['dyntemp']*kelvin, \
                                    options['md']['langfbeta']/picosecond, \
                                    options['md']['dyntstep']*picosecond)
    #
    simulation = Simulation(psf.topology, sys, integrator, platform, properties)
    simulation.context.setPositions(crd.positions)
    with open(chk_fn, 'rb') as fp:
        simulation.context.loadCheckpoint(fp.read())

    simulation.reporters.append(StateDataReporter('%s.equil.log'%output_prefix, 2500, step=True, \
        time=True, kineticEnergy=True, potentialEnergy=True, temperature=True, progress=True, \
        remainingTime=True, speed=True, volume=True, totalSteps=steps_left, separator='\t'))
        #
    simulation.step(steps_left)
    #
    state = simulation.context.getState(getPositions=True,\
                                        getVelocities=True,\
                                        getForces=True,\
                                        getEnergy=True, \
                                        enforcePeriodicBox=True)
    #
    chk_fn = '%s.equil.restart'%(output_prefix)
    with open("%s.pkl"%chk_fn, 'wb') as fout:
        pickle.dump(state, fout)
    #
    boxinfo = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(nanometer)
    #
    equil_pdb_fn = path.Path('%s.equil.pdb'%output_prefix)
    pdb.xyz = state.getPositions(asNumpy=True).value_in_unit(nanometer)[None,:]
    pdb.unitcell_vectors = boxinfo[None,:]
    pdb.save(equil_pdb_fn.short())
    #
    cmd = []
    cmd.append("%s/update_water_name.py"%EXEC_HOME)
    cmd.append(equil_pdb_fn.short())
    system(cmd, verbose=verbose)
    #
    if 'use_modified_CMAP' in options['ff'] and options['ff']['use_modified_CMAP']:
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
    if 'ligand' in options:
        init_pdb = path.Path("%s.init.pdb"%output_prefix)
        add_ligand(options['ligand'], input_pdb, init_pdb)
        update_ligand_name(init_pdb, options['ligand']['ligand_s'])
    else:
        init_pdb = input_pdb
    #
    pdb = mdtraj.load(init_pdb.short())
    update_residue_name(init_pdb, pdb)
    if options['md']['solvate'] > 0.0:
        orient_fn, solv_fn = solvate_pdb(output_prefix, pdb, options, verbose)
        tempfile_s.extend([orient_fn, solv_fn])
    else:
        solv_fn = init_pdb
    if 'ligand' in options:
        update_ligand_name(orient_fn, options['ligand']['ligand_s'])
        update_ligand_name(solv_fn, options['ligand']['ligand_s'])
    #
    psf_fn, crd_fn = generate_PSF(output_prefix, solv_fn, options, verbose)
    tempfile_s.append(crd_fn)

    if 'ligand' in options:
        ligand_restraint = get_ligand_restratint(pdb, psf_fn, options['ligand'])
        options['ff']['ligand'] = ligand_restraint
    #
    equil_md(output_prefix, solv_fn, psf_fn, crd_fn, options, verbose)
    if 'ligand' in options:
        update_ligand_name(solv_fn, options['ligand']['ligand_s'])

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
        for key in options:
            options[key] = JSONdeserialize(options[key])
    if 'ligand_json' in options:
        options['ligand'] = read_ligand_json(options['ligand_json'])
    #
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
