#!/usr/bin/env python

import os
import sys
import path
import argparse
import numpy as np
import subprocess as sp
from importlib import import_module

from libcommon import *
from libmain import *

N_MODEL = 5


def main():
    arg = argparse.ArgumentParser(prog="trRosetta+PREFMD")
    arg.add_argument(dest="title", help="Job title")
    arg.add_argument("-i", "--input", dest="input_fa", help="input FA file")
    arg.add_argument(
        "-d", "--dir", dest="work_dir", default="./", help="working directory (default=./)"
    )
    arg.add_argument(
        "--keep",
        dest="keep",
        action="store_true",
        default=False,
        help="set temporary file mode (default=False)",
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
    arg.add_argument(
        "--hybrid", dest="use_hybrid", action="store_true", default=False, help="use hybrid"
    )
    arg.add_argument("--membrane", dest="is_membrane_protein", action="store_true", default=False)
    arg.add_argument("--ligand", dest="has_ligand", action="store_true", default=False)
    arg.add_argument("--oligomer", dest="is_oligomer", action="store_true", default=False)

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    if arg.input_fa is not None:
        arg.input_fa = path.Path(arg.input_fa)
    #
    # init
    if arg.title.endswith("job.json"):
        input_json = path.Path(arg.title)
        job = Job.from_json(input_json)
        job.verbose = arg.verbose
        job.keep_tmp = arg.keep
        job.append_to_joblist()
    else:
        job = import_module("init_sp").prep(arg)
        job.run_type = "trRosetta"
        job.run_exec = os.path.abspath(__file__)

    # trRosetta
    import_module("trRosetta").prep(job, job.init_fa)
    if not run(job, arg.wait_after_run):
        return
    #
    job.remove_from_joblist()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit()
