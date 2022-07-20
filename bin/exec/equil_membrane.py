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

WORK_HOME = os.getenv("PIPE_HOME")
assert WORK_HOME is not None
sys.path.insert(0, "%s/bin" % WORK_HOME)

import path
from libcommon import *

from libcustom import *
from libmd import construct_restraint, construct_membrane_restraint


def equil_md(output_prefix, pdb, psf_fn, crd_fn, options, verbose):
    psf = CharmmPsfFile(psf_fn.short())
    crd = CharmmCrdFile(crd_fn.short())
    crd.positions += options["translate"] * nanometer
    #
    box = pdb.atom_slice(pdb.top.select("resname HOH")).xyz[0]
    # box = np.array(crd.positions.value_in_unit(nanometers), dtype=float)
    boxsize = np.max(box, 0) - np.min(box, 0) + 0.3
    psf.setBox(*boxsize)
    #
    ff = CharmmParameterSet(*options["ff"]["toppar"])
    platform = Platform.getPlatformByName(options["openmm"]["platform"])
    #
    sys = psf.createSystem(
        ff,
        nonbondedMethod=PME,
        switchDistance=0.8 * nanometers,
        nonbondedCutoff=1.0 * nanometers,
        constraints=HBonds,
    )
    #
    sys.addForce(construct_restraint(psf, pdb, 0.5))
    sys.addForce(construct_membrane_restraint(psf, pdb, 0.1))
    #
    if "custom" in options["ff"] and options["ff"]["custom"] is not None:
        custom_restrains = read_custom_restraint(options["ff"]["custom"])
        custom_s = construct_custom_restraint(pdb, custom_restraints[1])
        for custom in custom_s:
            sys.addForce(custom)
    #
    steps_left = options["md"]["equil"][0]
    steps_heat = options["md"]["heat"][0]
    temp = options["md"]["heat"][1]
    temp_incr = options["md"]["heat"][2]
    #
    i = 0
    while temp < options["md"]["dyntemp"]:
        integrator = LangevinIntegrator(
            temp * kelvin,
            options["md"]["langfbeta"] / picosecond,
            options["md"]["dyntstep"] * picosecond,
        )
        #
        simulation = Simulation(psf.topology, sys, integrator, platform)
        simulation.context.setPositions(crd.positions)
        if i == 0:
            # state = simulation.context.getState(getEnergy=True)
            # print (state.getPotentialEnergy())

            simulation.minimizeEnergy(maxIterations=500)
            # state = simulation.context.getState(getEnergy=True)
            # print (state.getPotentialEnergy())
            simulation.context.setVelocitiesToTemperature(temp * kelvin)
        else:
            with open(chk_fn, "rb") as fp:
                simulation.context.loadCheckpoint(fp.read())
        simulation.reporters.append(
            StateDataReporter(
                "%s.heat.%03d.log" % (output_prefix, temp),
                500,
                step=True,
                time=True,
                kineticEnergy=True,
                potentialEnergy=True,
                temperature=True,
                progress=True,
                remainingTime=True,
                speed=True,
                volume=True,
                totalSteps=steps_heat,
                separator="\t",
            )
        )
        #
        simulation.step(steps_heat)
        #
        chk_fn = "%s.heat.restart" % (output_prefix)
        with open(chk_fn, "wb") as fout:
            fout.write(simulation.context.createCheckpoint())
        simulation = None  # have to free CUDA-related variables
        temp += temp_incr
        i += 1
        steps_left -= steps_heat
    #
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
    integrator = LangevinIntegrator(
        options["md"]["dyntemp"] * kelvin,
        options["md"]["langfbeta"] / picosecond,
        options["md"]["dyntstep"] * picosecond,
    )
    #
    simulation = Simulation(psf.topology, sys, integrator, platform)
    simulation.context.setPositions(crd.positions)
    with open(chk_fn, "rb") as fp:
        simulation.context.loadCheckpoint(fp.read())

    simulation.reporters.append(
        StateDataReporter(
            "%s.equil.log" % output_prefix,
            2500,
            step=True,
            time=True,
            kineticEnergy=True,
            potentialEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            volume=True,
            totalSteps=steps_left,
            separator="\t",
        )
    )
    #
    simulation.step(steps_left)
    #
    state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getForces=True,
        getEnergy=True,
        enforcePeriodicBox=True,
    )
    #
    boxinfo = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(nanometer)
    #
    equil_pdb_fn = path.Path("%s.equil.pdb" % output_prefix)
    pdb.xyz = state.getPositions(asNumpy=True).value_in_unit(nanometer)[None, :]
    pdb.unitcell_vectors = boxinfo[None, :]
    pdb.save(equil_pdb_fn.short())
    #
    cmd = []
    cmd.append("%s/update_water_name.py" % EXEC_HOME)
    cmd.append(equil_pdb_fn.short())
    system(cmd, verbose=verbose)
    #
    if "use_modified_CMAP" in options["ff"] and options["ff"]["use_modified_CMAP"]:
        cmd = ["%s/resName_modified_CMAP.py" % EXEC_HOME, equil_pdb_fn.short()]
        output = system(cmd, verbose=verbose, stdout=True)
        with equil_pdb_fn.open("wt") as fout:
            if "ssbond" in options:
                for line in options["ssbond"]:
                    fout.write("%s\n" % line)
            fout.write(output)
    #
    chk_fn = "%s.equil.restart" % (output_prefix)
    with open("%s.pkl" % chk_fn, "wb") as fout:
        pickle.dump(state, fout)


def run(input_pdb, input_psf, input_crd, output_prefix, options, verbose, nonstd):
    tempfile_s = []
    #
    pdb = mdtraj.load(input_pdb.short())
    translate = -np.min(pdb.xyz[0], 0)
    pdb.xyz += translate
    options["translate"] = translate
    #
    psf_fn = path.Path("%s.psf" % output_prefix)
    if not psf_fn.status():
        system(["cp", input_psf.short(), psf_fn.short()])
    orient_fn = path.Path("%s.orient.pdb" % output_prefix)
    if not orient_fn.status():
        prot = pdb.atom_slice(pdb.top.select("protein"))
        prot.save(orient_fn.short())
    #
    equil_md(output_prefix, pdb, input_psf, input_crd, options, verbose)


def main():
    arg = argparse.ArgumentParser(prog="equil")
    arg.add_argument(dest="output_prefix")
    arg.add_argument(dest="input_pdb")
    arg.add_argument(dest="input_psf")
    arg.add_argument(dest="input_crd")
    arg.add_argument("--input", dest="input_json", required=True)
    arg.add_argument("--toppar", dest="toppar", nargs="*", default=None)
    arg.add_argument("--custom", dest="custom_file", default=None)
    arg.add_argument("--temp", dest="temp", default=None, type=float)
    arg.add_argument("--non_standard", dest="non_standard", default=False, action="store_true")
    arg.add_argument("-v", "--verbose", default=False, action="store_true")
    arg.add_argument(
        "--keep",
        dest="keep",
        action="store_true",
        default=False,
        help="set temporary file mode (default=False)",
    )
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    #
    input_pdb = path.Path(arg.input_pdb)
    input_psf = path.Path(arg.input_psf)
    input_crd = path.Path(arg.input_crd)
    #
    with open(arg.input_json) as fp:
        options = json.load(fp)
    if arg.toppar is not None:
        options["ff"]["toppar"] = arg.toppar
    if arg.custom_file is not None:
        options["ff"]["custom"] = arg.custom_file
    if arg.temp is not None:
        options["md"]["dyntemp"] = arg.temp
    #
    run(input_pdb, input_psf, input_crd, arg.output_prefix, options, arg.verbose, arg.non_standard)


if __name__ == "__main__":
    main()
