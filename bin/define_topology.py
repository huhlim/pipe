#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

def prep(job, input_pdb):
    job.top_fn = job.init_home.fn("solute.pdb")
    n_atom = 0
    if (not job.top_fn.status()) or (not job.has("n_atom") or not job.has("ssbond")):
        pdb = [] ; ssbond = []
        with input_pdb.open() as fp:
            for line in fp:
                if line.startswith("ATOM"):
                    n_atom += 1
                    pdb.append(line)
                elif line.startswith("SSBOND"):
                    pdb.append(line)
                    ssbond.append(line.rstrip())
        with job.top_fn.open("wt") as fout:
            fout.writelines(pdb)
            fout.write("TER\n")
            fout.write("END\n")

    job.n_atom = n_atom
    job.ssbond = ssbond
    #
    job.to_json()

def main():
    arg = argparse.ArgumentParser(prog='define_topology')
    arg.add_argument(dest='command', choices=['prep'], help='exec type')
    arg.add_argument(dest='work_dir', help='work_dir, which has a JSON file')
    arg.add_argument('-i', '--input', dest='input_pdb', nargs='?', required=True, \
            help='input PDB file')

    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args()
    #
    if arg.work_dir.endswith(".json"):
        arg.json_job = path.Path(arg.work_dir)
    else:
        arg.work_dir = path.Dir(arg.work_dir)
        arg.json_job = arg.work_dir.fn("job.json")
    #
    job = Job.from_json(arg.json_job)
    #
    if arg.command == 'prep':
        arg.input_pdb = path.Path(arg.input_pdb)
        prep(job, arg.input_pdb)

if __name__ == '__main__':
    main()
