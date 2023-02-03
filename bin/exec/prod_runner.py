#!/usr/bin/env python

import os
import sys
import time
import json
import argparse
import pickle

import numpy as np

import mdtraj

from openmm.unit import *
from openmm.openmm import *
from openmm.app import *
from openmm.app.internal.charmm.exceptions import CharmmPSFWarning
#
from importlib import import_module

#
import warnings

warnings.filterwarnings("ignore", category=CharmmPSFWarning)

from libcustom import *
from libmd import construct_ligand_restraint
from libligand import read_ligand_json, get_ligand_restratint

import warnings

warnings.filterwarnings("ignore")


def run(arg, options):
    t_init = time.time()
    #
    pdb = mdtraj.load(arg.init_pdb)
    crd = pdb.xyz[0]
    box = pdb.unitcell_lengths[0]

    psf = CharmmPsfFile(options["input"]["psf"])
    psf.setBox(*box)
    #
    ff_file_s = options["ff"]["toppar"]
    if "ligand" in options and len(options["ligand"]["str_fn_s"]) > 0:
        ff_file_s.extend(options["ff"]["cgenff"])
        ff_file_s.extend([fn.short() for fn in options["ligand"]["str_fn_s"]])

    ff = CharmmParameterSet(*ff_file_s)
    properties = {}
    cuda_devices = os.getenv("CUDA_VISIBLE_DEVICES")
    if cuda_devices is not None:
        n_gpu = len(cuda_devices.split(","))
    else:
        n_gpu = 0
        options["openmm"]["platform"] = "CPU"
    if n_gpu > 1:
        properties["DeviceIndex"] = ",".join(["%d" % index for index in range(n_gpu)])
    platform = Platform.getPlatformByName(options["openmm"]["platform"])
    if options["openmm"]["platform"] == "CPU":
        properties = {}
    #
    switchDistance = options["ff"].get("switchDistance", 0.8) * nanometers
    nonbondedCutoff = options["ff"].get("nonbondedCutoff", 1.0) * nanometers
    #
    if arg.use_hmr:
        sys = psf.createSystem(
            ff,
            nonbondedMethod=PME,
            switchDistance=switchDistance,
            nonbondedCutoff=nonbondedCutoff,
            hydrogenMass=3.008 * amu,
            constraints=HBonds,
        )
    else:
        sys = psf.createSystem(
            ff,
            nonbondedMethod=PME,
            switchDistance=switchDistance,
            nonbondedCutoff=nonbondedCutoff,
            constraints=HBonds,
        )

    if arg.rsr_fn is not None:
        ref_fn, custom_rsr = read_custom_restraint(arg.rsr_fn)
        ref = mdtraj.load(ref_fn)
        rsr_s = construct_custom_restraint(ref, custom_rsr)
        for rsr in rsr_s:
            sys.addForce(rsr)
    #
    if "ligand" in options:
        ligand_restraint = get_ligand_restratint(pdb, options["input"]["psf"], options["ligand"])
        ligand_restraints = construct_ligand_restraint(ligand_restraint)
        sys.addForce(ligand_restraints)
    #
    if arg.barostat is None:
        pass
    elif arg.barostat == "isotropic":
        sys.addForce(MonteCarloBarostat(1.0 * bar, options["md"]["dyntemp"] * kelvin))
    elif arg.barostat == "anisotropic":
        sys.addForce(
            MonteCarloAnisotropicBarostat(
                [1.0 * bar, 1.0 * bar, 1.0 * bar], options["md"]["dyntemp"] * kelvin
            )
        )
    elif arg.barostat == "XYisotropic":
        sys.addForce(
            MonteCarloMembraneBarostat(
                1.0 * bar,
                0.0 * bar * nanometers,
                options["md"]["dyntemp"] * kelvin,
                MonteCarloMembraneBarostat.XYIsotropic,
                MonteCarloMembraneBarostat.ZFree,
                100,
            )
        )
    else:
        sys.exit("ERROR: unknown barostat.\n")
    #
    if arg.integrator == "Langevin":
        integrator = LangevinIntegrator(
            options["md"]["dyntemp"] * kelvin,
            options["md"]["langfbeta"] / picosecond,
            options["md"]["dyntstep"] * picosecond,
        )
    elif arg.integrator == "Berendsen":
        from libmd import BerendsenVelocityVerletIntegrator
        integrator = BerendsenVelocityVerletIntegrator(
            options["md"]["dyntemp"] * kelvin, options["md"]["dyntstep"] * picosecond
        )
    else:
        raise ValueError(arg.integrator)
    #
    simulation = Simulation(psf.topology, sys, integrator, platform, properties)
    simulation.context.setPositions(crd)
    #
    if arg.restart_fn is not None:
        with open(arg.restart_fn, "rb") as fp:
            restart_state = pickle.load(fp)
        simulation.context.setState(restart_state)
    else:
        simulation.context.setVelocitiesToTemperature(options["md"]["dyntemp"] * kelvin)
    simulation.context.setTime(0.0)
    simulation.context.setStepCount(0)
    #
    if "n_atom" in options["input"]:
        simulation.reporters.append(
            mdtraj.reporters.DCDReporter(
                arg.out_dcd_fn,
                options["md"]["dynoutfrq"],
                atomSubset=np.arange(options["input"]["n_atom"]),
            )
        )
    else:
        simulation.reporters.append(DCDReporter(arg.out_dcd_fn, options["md"]["dynoutfrq"]))

    if arg.out_log_fn is not None:
        simulation.reporters.append(
            StateDataReporter(
                arg.out_log_fn,
                options["md"]["dynoutfrq"],
                step=True,
                time=True,
                kineticEnergy=True,
                potentialEnergy=True,
                temperature=True,
                progress=True,
                remainingTime=True,
                speed=True,
                volume=True,
                totalSteps=options["md"]["dynsteps"],
                separator="\t",
            )
        )
    #
    simulation.step(options["md"]["dynsteps"])
    #
    state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getForces=True,
        getEnergy=True,
        enforcePeriodicBox=True,
    )
    with open(arg.out_chk_fn, "wb") as fout:
        pickle.dump(state, fout)
    #
    t_final = time.time()
    t_spend = t_final - t_init
    #
    if arg.out_log_fn is not None:
        with open(arg.out_log_fn, "at") as fout:
            fout.write("ELAPSED TIME:       %8.2f SECONDS\n" % t_spend)
    else:
        sys.stdout.write("ELAPSED TIME:       %8.2f SECONDS\n" % t_spend)


def main():
    arg = argparse.ArgumentParser(prog="prod")
    arg.add_argument("--input", dest="input_json", required=True)
    arg.add_argument("--pdb", dest="init_pdb", required=True)
    arg.add_argument("--dcd", dest="out_dcd_fn", required=True)
    arg.add_argument("--chk", dest="out_chk_fn", required=True)
    arg.add_argument("--log", dest="out_log_fn")
    arg.add_argument("--custom", dest="rsr_fn")
    arg.add_argument("--restart", dest="restart_fn")
    arg.add_argument("--hmr", dest="use_hmr", default=False, action="store_true")
    arg.add_argument("--barostat", dest="barostat", default=None)
    arg.add_argument("--integrator", dest="integrator", default="Langevin")
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    with open(arg.input_json) as fp:
        options = json.load(fp)
    #
    if "ligand_json" in options:
        options["ligand"] = read_ligand_json(options["ligand_json"])

    run(arg, options)


if __name__ == "__main__":
    main()
