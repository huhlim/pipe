#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

def main():
    arg = argparse.ArgumentParser(prog='set_membrane_protein')
    arg.add_argument(dest='work_dir', help='work_dir, which has a JSON file')
    arg.add_argument('--pdb', dest='pdb_fn', nargs='*')
    arg.add_argument('--psf', dest='psf_fn', nargs='*')
    arg.add_argument('--crd', dest='crd_fn', nargs='*')

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
    job.is_membrane_protein = True
    job.membrane_pdb = [path.Path(fn) for fn in arg.pdb_fn]
    job.membrane_psf = [path.Path(fn) for fn in arg.psf_fn]
    job.membrane_crd = [path.Path(fn) for fn in arg.crd_fn]
    #
    job.to_json()

if __name__ == '__main__':
    main()
