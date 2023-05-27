#!/usr/bin/env python

import os
import sys
import path
import argparse
import subprocess as sp
from importlib import import_module

from libcommon import *
from libmain import *
from libligand import get_ligand_info


def main():
    arg = argparse.ArgumentParser(prog="PREFMD")
    arg.add_argument(dest="title", help="Job title")
    arg.add_argument("-i", "--input", dest="input_pdb", help="input PDB file")
    arg.add_argument(
        "-d", "--dir", dest="work_dir", default="./", help="working directory (default=./)"
    )
    arg.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="set verbose mode (default=False)",
    )
    arg.add_argument(
        "-w",
        "--wait",
        dest="wait_after_run",
        action="store_true",
        default=False,
        help="set running type (default=False)",
    )
    arg.add_argument("--membrane", dest="is_membrane_protein", action="store_true", default=False)
    arg.add_argument("--ligand", dest="has_ligand", action="store_true", default=False)
    arg.add_argument("--temperature", dest="temperature", default=298.15, type=float)
    arg.add_argument("--dynoutfrq", dest="dynoutfrq", default=0.1, type=float)
    arg.add_argument("--dynsteps", dest="dynsteps", default=20.0, type=float)
    arg.add_argument("--dyniter", dest="dyniter", default=10, type=int)
    arg.add_argument("--no-solute", dest="generate_solute_file", default=True, action="store_false")

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    if arg.input_pdb is not None:
        arg.input_pdb = path.Path(arg.input_pdb)

    # init
    if arg.title.endswith("job.json"):
        input_json = path.Path(arg.title)
        job = Job.from_json(input_json)
        job.verbose = arg.verbose
        job.keep_tmp = False
        job.append_to_joblist()
    else:
        job = import_module("init_md").prep(arg)
        job.run_type = "simulate_for_msm"
        job.run_exec = path.Path(__file__).path()
        #
        dt = 0.002
        job.dynoutfrq = int(arg.dynoutfrq / dt * 1000)
        job.dynsteps = int(arg.dynsteps / dt * 1000)
        job.temperature = arg.temperature
        job.dyniter = arg.dyniter
        job.generate_solute_file = arg.generate_solute_file
    job.to_json()
    #
    runner_s = [job.init_pdb[0]]
    for fn in job.init_home.glob("*.pdb"):
        if fn in runner_s:
            continue
        if fn.name() in ["init", "solute", "segname", "complete"]:
            continue
        runner_s.append(fn)
    n_init = len(runner_s)

    # define topology
    import_module("define_topology").prep(job, job.init_pdb[0])
    if job.has("is_membrane_protein"):
        if not get_membrane_topology(job, n_init, arg.wait_after_run):
            sys.stderr.write("waiting for CHARMM-GUI membrane topology... \n")
            return
    if job.has("has_ligand"):
        if not get_ligand_info(job, arg.wait_after_run):
            sys.stderr.write("waiting for ligand info... \n")
            return
    job.to_json()
    #
    json_update = {}
    json_update["md"] = {
        "dyntemp": job.temperature,
    }

    # equil
    if not job.has("is_membrane_protein"):
        import_module("equil").prep(
            job,
            0,
            runner_s,
            path.Path("%s/equil_for_msm.json" % (DEFAULT_HOME), update=json_update),
        )
    else:
        import_module("equil").prep_membrane(
            job,
            0,
            job.membrane_pdb,
            job.membrane_psf,
            job.membrane_crd,
            path.Path("%s/equil_membrane_for_msm.json" % (DEFAULT_HOME)),
            update=json_update,
        )
    if not run(job, arg.wait_after_run):
        return

    json_update["md"].update(
        {
            "dynoutfrq": job.dynoutfrq,
            "dynsteps": job.dynsteps,
            "iter": job.dyniter,
        }
    )
    if not arg.generate_solute_file:
        json_update["generate_solute_file"] = False

    # prod
    n_traj_per_init = 10
    if not job.has("is_membrane_protein"):
        prod_input = path.Path("%s/prod_for_msm.json" % DEFAULT_HOME)
        for i in range(n_init):
            import_module("prod").prep(job, i, i, prod_input, n_traj_per_init, update=json_update)
    else:
        prod_input = path.Path("%s/prod_membrane_for_msm.json" % DEFAULT_HOME)
        for i in range(n_init):
            import_module("prod").prep(job, i, i, prod_input, n_traj_per_init, update=json_update)
    if not run(job, arg.wait_after_run):
        return


if __name__ == "__main__":
    main()
