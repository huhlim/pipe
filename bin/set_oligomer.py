#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *


def main():
    arg = argparse.ArgumentParser(prog="set_oligomer")
    arg.add_argument(dest="work_dir", help="work_dir, which has a JSON file")
    arg.add_argument("--unset", dest="unset_oligomer", default=False, action="store_true")

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
    if arg.unset_oligomer:
        if job.has("is_oligomer"):
            del job.is_oligomer
    else:
        job.is_oligomer = True
    #
    job.to_json()


if __name__ == "__main__":
    main()
