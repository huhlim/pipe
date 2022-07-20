#!/usr/bin/env python

import os
import sys
import path
import json
import argparse
import numpy as np

from libcommon import *


def prep(arg):
    work_home = path.Dir("%s/%s" % (arg.work_dir, arg.title))
    json_job = work_home.fn("job.json")
    if json_job.status():
        job = Job.from_json(json_job)
        job.verbose = arg.verbose
        job.keep_tmp = arg.keep
        job.append_to_joblist()
        return job
    #
    assert arg.input_pdb is not None
    #
    job = Job(arg.work_dir, arg.title, build=True)
    job.run_type = "refine_meta"
    #
    job.init_home = job.work_home.subdir("init", build=True)
    job.verbose = arg.verbose
    job.keep_tmp = arg.keep
    #
    out = job.init_home.fn("init.pdb")
    if not out.status():
        cmd = ["convpdb.pl", "-out", "generic", arg.input_pdb.short()]
        output = system(cmd, stdout=True, verbose=job.verbose)
        with out.open("wt") as fout:
            fout.write(output)
    job.init_pdb = [out]
    job.to_json()
    job.append_to_joblist()
    return job
