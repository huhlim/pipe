#!/usr/bin/env python

import os
import sys
import path
import argparse
import subprocess as sp
from importlib import import_module

from libcommon import *
from libmain import *

EXEC_REFINE = '%s/casp14_refine.py'%BIN_HOME
N_MODEL_REFINE = 5

def run_refine(title, input_pdb, work_home, **kwargs):
    cmd = []
    cmd.append(EXEC_REFINE)
    cmd.append(title)
    cmd.extend(['--input', input_pdb.short()])
    cmd.extend(['--dir', work_home.short()])
    if kwargs.get("verbose", False):
        cmd.append("--verbose")
    if kwargs.get("wait_after_run", False):
        cmd.append("--wait")
    if kwargs.get("use_hybrid", False):
        cmd.append("--hybrid")
    #
    proc = sp.Popen(cmd)
    return proc

def wait_refine(refine_proc_s, refine_home_s, wait_after_run, sleep=30):
    # check proc_s
    if refine_proc_s is not None:
        while True:
            proc_status = True
            for refine_proc in refine_proc_s:
                if refine_proc.poll() is None:
                    proc_status = False
            #
            if proc_status:
                break
            #
            time.sleep(sleep)
    #
    # check model_s
    status = True
    refined = []
    for refine_home in refine_home_s:
        model_home = refine_home.subdir("final")
        model_s = model_home.glob("model_?.pdb")
        refined.append(model_s)
        if len(model_s) < N_MODEL_REFINE:
            status = False
    #
    if wait_after_run and status:
        # expect to have refined models after all the process ends
        # , so if not raise error.
        sys.exit("Error: failed to get refined models\n")

    return status, refined

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
    if arg.title.endswith("/job.json"):
        input_json = path.Path(arg.title)
        job = Job.from_json(input_json)
        job.verbose = arg.verbose
        job.keep_tmp = arg.keep
        job.append_to_joblist()
    else:
        job = import_module("init_sp").prep(arg)
    
    # trRosetta
    import_module("trRosetta").prep(job, job.init_fa)
    if not run(job, arg.wait_after_run):
        return 
    domain_pdb_s = get_outputs(job, 'trRosetta', expand='model_s')#

    # create refine directory
    job.refine_home = job.work_home.subdir("refine", build=True)
    #
    if not job.has("refine_s"): # first time 
        refine_proc_s = []
        job.refine_s = []
        for out in domain_pdb_s:
            pdb_fn = out[0]
            domain_id = pdb_fn.name()
            #
            refine_proc = run_refine(domain_id, pdb_fn, job.refine_home, verbose=arg.verbose, \
                                     use_hybrid=arg.use_hybrid, wait_after_run=arg.wait_after_run)
            refine_home = job.refine_home.subdir(domain_id)
            #
            refine_proc_s.append(refine_proc)
            job.refine_s.append(refine_home)
    else:
        refine_proc_s = None
    job.to_json()
    #
    status, refined = wait_refine(refine_proc_s, job.refine_s, arg.wait_after_run)
    if not status:
        return
    #
    # TODO: paste refined models
    #
    # locPREFMD
    #
    # final

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit()

