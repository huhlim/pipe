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
    job = Job(arg.work_dir, arg.title, build=True)
    job.run_type = 'meta'
    job.init_home = job.work_home.subdir("init", build=True)
    job.meta_json_s = arg.meta_json_s
    job.verbose = arg.verbose
    job.keep_tmp = arg.keep
    #
    job.to_json()
    job.append_to_joblist()
    return job

def override(arg):
    pass

def main():
    arg = argparse.ArgumentParser(prog='init')
    arg.add_argument(dest='command', choices=['prep', 'override'], help='exec type')
    arg.add_argument(dest='title', help='Job title')
    arg.add_argument('-i', '--input', dest='meta_json_s', required=True, nargs='*', \
            help='input JSON files')
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
        arg.meta_json_s = [path.Path(json_fn) for json_fn in arg.meta_json_s]
        prep(arg)
    elif arg.command == 'override':
        override(arg)

if __name__ == '__main__':
    main()
