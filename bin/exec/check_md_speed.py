#!/usr/bin/env python

import os
import sys
import json
import time
import argparse
import pickle
from tempfile import TemporaryDirectory

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import scipy.spatial.distance

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

def check_speed(output_prefix, solv_fn, psf_fn, crd_fn, restart_fn, options, verbose):
    psf = CharmmPsfFile(psf_fn.short())
    pdb = mdtraj.load(solv_fn.short())
    #
    if crd_fn is not None:
        crd = CharmmCrdFile(crd_fn.short())
        crd = crd.positions.value_in_unit(nanometers)
        box = np.array(crd, dtype=float)
        boxsize = np.max(box, 0) - np.min(box, 0)
    else:
        crd = pdb.xyz[0]
        boxsize = pdb.unitcell_lengths[0]
    psf.setBox(*boxsize)
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
    steps = 25000
    temp = 10
    #
    integrator = LangevinIntegrator(temp*kelvin, 1.0/picosecond, 0.001*picosecond)
    #
    simulation = Simulation(psf.topology, sys, integrator, platform)
    simulation.context.setPositions(crd)
    #
    if restart_fn is None:
        simulation.minimizeEnergy(maxIterations=500)
        simulation.context.setVelocitiesToTemperature(temp*kelvin)
    else:
        with open(restart_fn.short(), 'rb') as fp:
            restart_state = pickle.load(fp)
        simulation.context.setState(restart_state)
    #
    t0 = time.time()
    simulation.step(steps)
    t1 = time.time()
    #
    t_spend = t1-t0
    speed = float(steps) / t_spend  # steps / s
    return speed

def run(output_prefix, options, verbose, nonstd):
    cwd = os.getcwd()
    tmpdir = TemporaryDirectory(prefix='md_speed.')
    os.chdir(tmpdir.name)
    #
    if 'input' in options:  # prod
        solv_fn = options['input']['pdb']
        psf_fn = options['input']['psf']
        crd_fn = None
        restart_fn = options['input']['restart']
    else:
        pdb = mdtraj.load(options['input_pdb'].short())
        orient_fn, solv_fn = solvate_pdb(output_prefix, pdb, options, verbose)
        #
        psf_fn, crd_fn = generate_PSF(output_prefix, solv_fn, options, verbose)
        restart_fn = None
    #
    speed = check_speed(output_prefix, solv_fn, psf_fn, crd_fn, restart_fn, options, verbose)
    #
    os.chdir(cwd)
    with open("SPEED", 'wt') as fout:
        fout.write("%10.4f\n"%speed)

def main():
    arg = argparse.ArgumentParser(prog='check_md_speed')
    arg.add_argument(dest='output_prefix')
    arg.add_argument('--input', dest='input_json', required=True)
    arg.add_argument('--non_standard', dest='non_standard', default=False, action='store_true')
    arg.add_argument('-v', '--verbose', default=False, action='store_true')
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    #
    with open(arg.input_json) as fp:
        options = json.load(fp)
        if 'input' in options:
            for keyword in ['pdb', 'psf', 'restart']:
                options['input'][keyword] = path.Path(options['input'][keyword])
        elif 'input_pdb' in options:
            options['input_pdb'] = path.Path(options['input_pdb'])
    #
    run(arg.output_prefix, options, arg.verbose, arg.non_standard)

if __name__=='__main__':
    main()
