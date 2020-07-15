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
from casp14_sp import wait_refine, paste_refined, update_bfactor

EXEC_REFINE = '%s/casp14_refine_meta.py'%BIN_HOME

N_MODEL_REFINE = 5
N_MODEL = 5

def check_meta_tasks(meta_s, wait_after_run, sleep=30):
    while True:
        status = True
        #
        for meta in meta_s:
            method = 'prod'
            if len(meta.get_task(method)) == 0:
                status = False ; break

            method = 'prod'
            if len(meta.get_task(method, not_status='DONE')) > 0:
                status = False ; break
            #
            method = 'score'
            if len(meta.get_task(method)) == 0:
                status = False ; break
            #
            method = 'score'
            if len(meta.get_task(method, not_status='DONE')) > 0:
                status = False ; break
        #
        if status: break
        #
        if wait_after_run:
            time.sleep(sleep)
        else:
            break
    return status

def run_refine(title, input_pdb, work_home, sp_refine_job, refine_job_s, **kwargs):
    cmd = []
    cmd.append(EXEC_REFINE)
    cmd.append(title)
    cmd.extend(['--input', input_pdb.short()])
    cmd.extend(['--dir', work_home.short()])
    cmd.extend(['--sp', sp_refine_job.json_job.short()])
    cmd.append('--meta')
    cmd.extend([refine_job.json_job.short() for refine_job in refine_job_s])
    if kwargs.get("verbose", False):
        cmd.append("--verbose")
    if kwargs.get("wait_after_run", False):
        cmd.append("--wait")
    #
    proc = sp.Popen(cmd)
    time.sleep(60)
    return proc

def main():
    arg = argparse.ArgumentParser(prog='casp14_meta')
    arg.add_argument(dest='title', help='Job title')
    arg.add_argument('-i', '--input', dest='meta_json_s', nargs='*', \
            help='input JSON files', default=[])
    arg.add_argument('-d', '--dir', dest='work_dir', default='./',\
            help='working directory (default=./)')
    arg.add_argument('-w', '--wait', dest='wait_after_run', action='store_true', default=False,\
            help='set running type (default=False)')
    arg.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,\
            help='set verbose mode (default=False)')
    arg.add_argument('--keep', dest='keep', action='store_true', default=False,\
            help='set temporary file mode (default=False)')

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    arg.meta_json_s = [path.Path(fn) for fn in arg.meta_json_s]

    # init
    if arg.title.endswith("job.json"):
        input_json = path.Path(arg.title)
        job = Job.from_json(input_json)
        job.verbose = arg.verbose
        job.keep_tmp = arg.keep
        job.append_to_joblist()
    else:
        if len(arg.meta_json_s) == 0:
            sys.exit("Input json files were not given.\n")
        job = import_module("init_meta").prep(arg)
    #
    meta_s = []
    for meta_json_fn in job.meta_json_s:
        meta_s.append(Job.from_json(meta_json_fn))
    n_meta = len(meta_s)
    #
    sp_job = meta_s[0]
    sp_job.sub_s = []
    for refine_home in sp_job.refine_s:
        sp_job.sub_s.append(Job.from_json(refine_home.fn("job.json")))
    refine_job_s = meta_s[1:]
    #
    if not check_meta_tasks(sp_job.sub_s, arg.wait_after_run):
        return
    if not check_meta_tasks(refine_job_s, arg.wait_after_run):
        return
    #
    domain_pdb_s, trRosetta_min = get_outputs(sp_job, 'trRosetta', expand='model_s')[0]
    #
    # create refine directory
    job.refine_home = job.work_home.subdir("refine", build=True)
    #
    has_refine = job.has("refine_s")
    if has_refine:
        for refine_home in job.refine_s:
            if not refine_home.fn("job.json").status():
                has_refine = False ; break
    #
    if not has_refine:
        refine_proc_s = []
        job.refine_s = []
        for domain_pdb, sp_refine_job in zip(domain_pdb_s, sp_job.sub_s):
            domain_id = sp_refine_job.title
            #
            refine_proc = run_refine(domain_id, domain_pdb, job.refine_home, \
                                     sp_refine_job, refine_job_s, \
                                     verbose=job.verbose, wait_after_run=arg.wait_after_run)
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
    # paste
    job.work_home.chdir()
    model_home = job.work_home.subdir("model", build=True)
    prep_s = paste_refined(model_home, trRosetta_min, refined, verbose=job.verbose)
    #
    import_module("scwrl").prep(job, prep_s)
    if not run(job, arg.wait_after_run):
        return
    scwrl_out = get_outputs(job, 'scwrl')
    #
    import_module("locPREFMD").prep(job, [out[0] for out in scwrl_out])
    if not run(job, arg.wait_after_run):
        return 
    locPREFMD_out = get_outputs(job, "locPREFMD")

    # final
    final_home = job.work_home.subdir("final", build=True)
    for i,out in enumerate(locPREFMD_out):
        pdb_fn = final_home.fn("model_%d.pdb"%(i+1))
        if not pdb_fn.status():
            update_bfactor(pdb_fn, out[0], prep_s[i])
    #
    job.remove_from_joblist()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit()

