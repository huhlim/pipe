#!/usr/bin/env python

import os
import sys
import path
import json
import argparse
import numpy as np

from libcommon import *

FF = {"c36m": [f"{BIN_HOME}/ff/c36m/par_all36m_prot.prm", f"{BIN_HOME}/ff/c36m/top_all36_prot.rtf"]}


def prep(arg, ff=FF["c36m"]):
    work_home = path.Dir("%s/%s" % (arg.work_dir, arg.title))
    json_job = work_home.fn("job.json")
    if json_job.status():
        job = Job.from_json(json_job)
        job.verbose = arg.verbose
        job.keep_tmp = False
        if hasattr(arg, "water_pdb"):
            job.water_pdb = path.Path(arg.water_pdb)
        job.append_to_joblist()
        return job
    #
    assert arg.input_pdb is not None
    #
    job = Job(arg.work_dir, arg.title, build=True)
    job.run_type = "run_md"
    #
    job.init_home = job.work_home.subdir("init", build=True)
    job.verbose = arg.verbose
    job.keep_tmp = False
    if hasattr(arg, "water_pdb"):
        job.water_pdb = path.Path(arg.water_pdb)
    #
    seg_pdb = job.init_home.fn("segname.pdb")
    if not seg_pdb.status():
        cmd = [f"{EXEC_HOME}/put_segnames.py"]
        cmd.append(arg.input_pdb.short())
        cmd.append(seg_pdb.short())
        system(cmd, verbose=job.verbose)
    #
    complete_pdb = job.init_home.fn("complete.pdb")
    if not complete_pdb.status():
        cmd = [f"{EXEC_HOME}/complete.py"]
        cmd.append(seg_pdb.short())
        cmd.append("--toppar")
        cmd.extend(ff)
        # cmd.append("--scwrl")
        cmd.append("-o")
        cmd.append(complete_pdb.short())
        system(cmd, verbose=job.verbose)
    #
    init_pdb = job.init_home.fn("init.pdb")
    if not init_pdb.status():
        cmd = ["convert_to_generic.py", complete_pdb.short()]
        # cmd = ["convpdb.pl", "-out", "generic", complete_pdb.short()]
        output = system(cmd, stdout=True, verbose=job.verbose)
        with init_pdb.open("wt") as fout:
            fout.write(output)
    #
    job.init_pdb = [init_pdb]
    if arg.is_membrane_protein:
        job.is_membrane_protein = True
    if arg.has_ligand:
        job.has_ligand = True
    job.to_json()
    job.append_to_joblist()
    return job


def override(arg):
    pass


def main():
    arg = argparse.ArgumentParser(prog="init")
    arg.add_argument(dest="command", choices=["prep", "override"], help="exec type")
    arg.add_argument(dest="title", help="Job title")
    arg.add_argument("-i", "--input", dest="input_pdb", required=True, help="input PDB file")
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

    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args()
    #
    if arg.command == "prep":
        arg.input_pdb = path.Path(arg.input_pdb)
        prep(arg)
    elif arg.command == "override":
        override(arg)


if __name__ == "__main__":
    main()
