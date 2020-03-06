#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

METHOD = ''
EXEC = ''

def prep(job):
    pass

def run(job):
    pass

def submit(job):
    pass

def status(job):
    pass

def main():
    arg = argparse.ArgumentParser(prog='locPREFMD')
    arg.add_argument(dest='command', choices=['prep', 'run'], help='exec type')
    arg.add_argument(dest='work_dir', help='work_dir, which has a JSON file')
    arg.add_argument('-i', '--input', dest='input_pdb', nargs='?', \
            help='input PDB file, mandatory for "prep"')  

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
        if arg.input_pdb is None:
            sys.exit("Error: input_pdb required\n")
        arg.input_pdb = path.Path(arg.input_pdb)
        #
        prep(job, arg.input_pdb)

    elif arg.command == 'run':
        run(job)

if __name__ == '__main__':
    main()
