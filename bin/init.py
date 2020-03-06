#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

def prep(arg):
    job = Job(arg.work_dir, arg.title, build=True)
    job.init_home = job.work_home.subdir("init", build=True)
    job.verbose = arg.verbose
    job.keep_tmp = arg.keep
    #
    out = job.init_home.fn("init.pdb")
    if not out.status():
        cmd = ['convpdb.pl', '-out', 'generic', arg.input_pdb.short()]
        output = system(cmd, stdout=True, verbose=job.verbose)
        with out.open("wt") as fout:
            fout.write(output)
    job.init_pdb = [out]
    job.to_json()

def override(arg):
    pass

def main():
    arg = argparse.ArgumentParser(prog='init')
    arg.add_argument(dest='command', choices=['prep', 'override'], help='exec type')
    arg.add_argument(dest='title', help='Job title')
    arg.add_argument('-i', '--input', dest='input_pdb', nargs='?', required=True, \
            help='input PDB file')
    arg.add_argument('-d', '--dir', dest='work_dir', default='./',\
            help='working directory (default=./)')
    #arg.add_argument('-p', '--protocol', dest='protocol', default='standard', choices=['standard'], \
    #        help='protocol type (default=standard)')
    #arg.add_argument('--resource', dest='resource', default=None, \
    #        help='set resource via a JSON file')
    arg.add_argument('--keep', dest='keep', action='store_true', default=False,\
            help='set temporary file mode (default=False)')
    arg.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,\
            help='set verbose mode (default=False)')

    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args()
    #
    if arg.command == 'prep':
        arg.input_pdb = path.Path(arg.input_pdb)
        prep(arg)
    elif arg.command == 'override':
        override(arg)

if __name__ == '__main__':
    main()
