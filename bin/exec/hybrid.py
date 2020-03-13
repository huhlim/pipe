#!/usr/bin/env python

import os
import sys
import json
import argparse
from importlib import import_module

import warnings
warnings.filterwarnings("ignore")

WORK_HOME = os.getenv("PREFMD_HOME")
assert WORK_HOME is not None
sys.path.insert(0, '%s/bin'%WORK_HOME)

import path
from libcommon import *
from libpdb import Sequence, PDB

def run(arg):
    model_list = path.Path("hybrid/init_s")
    if not model_list.status():
        cmd = ['%s/prep_hybrid.py'%EXEC_HOME, arg.output_prefix, arg.init_pdb]
        system(cmd)
    #
    cmd = ['%s/run_hybrid.py'%EXEC_HOME, model_list.short()]
    system(cmd)

def main():
    arg = argparse.ArgumentParser(prog='hybrid')
    arg.add_argument(dest='output_prefix')
    arg.add_argument(dest='init_pdb')
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    #
    run(arg)

if __name__=='__main__':
    main()
