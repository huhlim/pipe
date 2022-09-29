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

from openmm.unit import *
from openmm.openmm import *
from openmm.app import *

WORK_HOME = os.getenv("PIPE_HOME")
assert WORK_HOME is not None
sys.path.insert(0, "%s/bin" % WORK_HOME)

import path
from libcommon import *
from libligand import read_ligand_json, add_ligand, update_ligand_name, get_ligand_restratint

from libcustom import *
from libmd import (
    solvate_pdb,
    generate_PSF,
    construct_restraint,
    construct_water_restraint,
    construct_ligand_restraint,
    update_residue_name,
)


def minimize(output_prefix, solv_fn, psf_fn, crd_fn, options, verbose):
    psf = CharmmPsfFile(psf_fn.short())
    crd = CharmmCrdFile(crd_fn.short())
    #
    pdb = mdtraj.load(solv_fn.short(), standard_names=False)
    #
    if options["md"].get("solvate", -1) > 0.0 and options["md"].get(
        "orient", True
    ):  # solvated via libmd.solvate_pdb
        box = np.array(crd.positions.value_in_unit(nanometers), dtype=float)
        boxsize = np.max(box, 0) - np.min(box, 0)
    else:  # solvated prior to the script
        boxsize = pdb.unitcell_lengths[0]
    psf.setBox(*boxsize)
    #
    ff_file_s = options["ff"]["toppar"]
    if "ligand" in options and len(options["ligand"]["str_fn_s"]) > 0:
        ff_file_s.extend(options["ff"]["cgenff"])
        ff_file_s.extend([fn.short() for fn in options["ligand"]["str_fn_s"]])

    ff = CharmmParameterSet(*ff_file_s)
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
    if options["md"].get("solvate_cryst", None) is not None:
        n_cryst_water = options["md"].get("n_cryst_water", -1)
        if n_cryst_water < 0:
            raise NotImplementedError
        sys.addForce(construct_water_restraint(psf, pdb, n_cryst_water, 0.5))
    #
    if "custom" in options["ff"] and options["ff"]["custom"] is not None:
        custom_restraints = read_custom_restraint(options["ff"]["custom"])
        custom_s = construct_custom_restraint(pdb, custom_restraints[1])
        for custom in custom_s:
            sys.addForce(custom)
    #
    if "ligand" in options["ff"]:
        ligand_restraints = construct_ligand_restraint(options["ff"]["ligand"])
        sys.addForce(ligand_restraints)
    #
    integrator = LangevinIntegrator(298.15 * kelvin, 1.0 / picosecond, 0.002 * picosecond)
    #
    simulation = Simulation(psf.topology, sys, integrator, platform)
    simulation.context.setPositions(crd.positions)

    state = simulation.context.getState(getEnergy=True)
    print(state.getPotentialEnergy())

    simulation.minimizeEnergy(maxIterations=500)

    state = simulation.context.getState(getEnergy=True)
    print(state.getPotentialEnergy())
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
    out_pdb_fn = path.Path("%s.pdb" % output_prefix)
    xyz = state.getPositions(asNumpy=True).value_in_unit(nanometer)[None, :]
    if pdb.xyz.shape[1] == xyz.shape[1]:
        pdb.xyz = xyz
        pdb.unitcell_vectors = boxinfo[None, :]
        pdb.save(out_pdb_fn.short())
    else:
        with out_pdb_fn.open("wt") as fout:
            PDBFile.writeFile(simulation.topology, state.getPositions(), fout)

    #
    cmd = []
    cmd.append("%s/update_water_name.py" % EXEC_HOME)
    cmd.append(out_pdb_fn.short())
    system(cmd, verbose=verbose)
    #
    if "use_modified_CMAP" in options["ff"] and options["ff"]["use_modified_CMAP"]:
        cmd = ["%s/resName_modified_CMAP.py" % EXEC_HOME, out_pdb_fn.short()]
        output = system(cmd, verbose=verbose, stdout=True)
        with out_pdb_fn.open("wt") as fout:
            if "ssbond" in options:
                for line in options["ssbond"]:
                    fout.write("%s\n" % line)
            fout.write(output)


def run(input_pdb, output_prefix, options, verbose):
    tempfile_s = []
    #
    if "ligand" in options:
        init_pdb = path.Path("%s.init.pdb" % output_prefix)
        add_ligand(options["ligand"], input_pdb, init_pdb)
        update_ligand_name(init_pdb, options["ligand"]["ligand_s"])
    else:
        init_pdb = input_pdb
    #
    pdb = mdtraj.load(init_pdb.short())
    update_residue_name(init_pdb, pdb)
    if options["md"].get("solvate", -1) > 0.0:
        orient_fn, solv_fn = solvate_pdb(output_prefix, pdb, options, verbose)
        tempfile_s.extend([orient_fn, solv_fn])
    else:
        solv_fn = init_pdb
    if "ligand" in options:
        update_ligand_name(orient_fn, options["ligand"]["ligand_s"])
        update_ligand_name(solv_fn, options["ligand"]["ligand_s"])
    #
    psf_fn, crd_fn = generate_PSF(output_prefix, solv_fn, options, verbose)
    tempfile_s.append(crd_fn)

    if "ligand" in options:
        ligand_restraint = get_ligand_restratint(pdb, psf_fn, options["ligand"])
        options["ff"]["ligand"] = ligand_restraint
    #
    minimize(output_prefix, solv_fn, psf_fn, crd_fn, options, verbose)
    if "ligand" in options:
        update_ligand_name(solv_fn, options["ligand"]["ligand_s"])


def main():
    arg = argparse.ArgumentParser(prog="equil")
    arg.add_argument(dest="output_prefix")
    arg.add_argument(dest="input_pdb")
    arg.add_argument("--input", dest="input_json", required=True)
    arg.add_argument("--toppar", dest="toppar", nargs="*", default=None)
    arg.add_argument("--custom", dest="custom_file", default=None)
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
    #
    with open(arg.input_json) as fp:
        options = json.load(fp)
        for key in options:
            options[key] = JSONdeserialize(options[key])
    if "ligand_json" in options:
        options["ligand"] = read_ligand_json(options["ligand_json"])
    #
    if arg.toppar is not None:
        options["ff"]["toppar"] = arg.toppar
    if arg.custom_file is not None:
        options["ff"]["custom"] = arg.custom_file
    #
    run(input_pdb, arg.output_prefix, options, arg.verbose)


if __name__ == "__main__":
    main()
