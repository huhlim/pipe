#!/usr/bin/env python

import os
import sys
import json
import time
import argparse
import subprocess as sp

import path
from libcommon import *

EXEC = {}
EXEC['refine'] = '%s/casp14_refine.py'%BIN_HOME
EXEC['sp'] = '%s/casp14_sp.py'%BIN_HOME

def run(verbose):
    with open("%s/bin/hosts/job_s.json"%WORK_HOME) as fp:
        job_fn_s = json.load(fp)
    #
    sys.stdout.write("TIME: %s\n"%time.ctime())
    #
    proc_s = []
    for job_fn in job_fn_s:
        job_fn = path.Path(job_fn)
        job = Job.from_json(job_fn)
        #
        cmd = [EXEC[job.run_type]]
        cmd.append(job_fn.path())
        #
        sys.stdout.write("PROC: %s\n"%(" ".join(cmd)))
        system(cmd, verbose=False)

def main():
    arg = argparse.ArgumentParser(prog='master_runner')
    arg.add_argument('--interval', dest='time_interval', default=60, type=int)
    arg.add_argument('--verbose', dest='verbose', action='store_true', default=False)
    #
    arg = arg.parse_args()
    #
    while True:
        run(arg.verbose)
        time.sleep(arg.time_interval)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit()
