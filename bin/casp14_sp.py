#!/usr/bin/env python

import os
import sys
import path
import argparse
import subprocess as sp
from importlib import import_module

from libcommon import *
from libmain import *

def main():
    arg = argparse.ArgumentParser(prog='trRosetta+PREFMD')
    arg.add_argument(dest='title', help='Job title')
    arg.add_argument('-i', '--input', dest='input_fa', \
            help='input FA file')
    arg.add_argument('-d', '--dir', dest='work_dir', default='./',\
            help='working directory (default=./)')
    arg.add_argument('--keep', dest='keep', action='store_true', default=False,\
            help='set temporary file mode (default=False)')
    arg.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,\
            help='set verbose mode (default=False)')
    arg.add_argument('-w', '--wait', dest='wait_after_run', action='store_true', default=False,\
            help='set running type (default=False)')
    arg.add_argument('--hybrid', dest='use_hybrid', action='store_true', default=False, \
            help='use hybrid')

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    if arg.input_fa is not None:
        arg.input_fa = path.Path(arg.input_fa)
    #
    # init
    job = import_module("init_sp").prep(arg)
    
    # trRosetta
    import_module("trRosetta").prep(job, job.init_fa)
    if not run(job, arg.wait_after_run):
        return 

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit()

