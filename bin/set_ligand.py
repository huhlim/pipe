#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *


def main():
    arg = argparse.ArgumentParser(prog="set_ligand")
    arg.add_argument(dest="work_dir", help="work_dir, which has a JSON file")
    arg.add_argument("--unset", dest="unset_ligand", default=False, action="store_true")

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
    if arg.unset_ligand:
        if job.has("has_ligand"):
            del job.has_ligand
    else:
        job.has_ligand = True
    #
    job.to_json()


if __name__ == "__main__":
    main()
