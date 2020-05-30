#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

def prep(arg):
    work_home = path.Dir("%s/%s"%(arg.work_dir, arg.title))
    json_job = work_home.fn("job.json")
    if json_job.status():
        job = Job.from_json(json_job)
        job.verbose = arg.verbose
        job.keep_tmp = arg.keep
        job.append_to_joblist()
        return job
    #
    assert arg.input_fa is not None
    #
    job = Job(arg.work_dir, arg.title, build=True)
    job.run_type = 'sp'
    job.use_hybrid = arg.use_hybrid
    job.verbose = arg.verbose
    job.keep_tmp = arg.keep
    if arg.is_membrane_protein:
        job.is_membrane_protein = True
    if arg.is_oligomer:
        job.is_oligomer = True
    #
    out = job.work_home.fn("%s.fa"%job.title)
    if not out.status():
        cmd = ['cp', arg.input_fa.short(), out.short()]
        system(cmd, verbose=job.verbose)
    job.init_fa = out
    job.to_json()
    job.append_to_joblist()
    return job

def override(arg):
    pass

def main():
    arg = argparse.ArgumentParser(prog='init')
    arg.add_argument(dest='command', choices=['prep', 'override'], help='exec type')
    arg.add_argument(dest='title', help='Job title')
    arg.add_argument('-i', '--input', dest='input_fa', required=True, \
            help='input PDB file')
    arg.add_argument('-d', '--dir', dest='work_dir', default='./',\
            help='working directory (default=./)')
    arg.add_argument('--keep', dest='keep', action='store_true', default=False,\
            help='set temporary file mode (default=False)')
    arg.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,\
            help='set verbose mode (default=False)')

    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args()
    #
    if arg.command == 'prep':
        arg.input_fa = path.Path(arg.input_fa)
        prep(arg)
    elif arg.command == 'override':
        override(arg)

if __name__ == '__main__':
    main()
