#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

METHOD = 'score'
EXEC = '%s/calc_statpot'%EXEC_HOME

def prep(job, input_dcd):
    if len(job.get_task(METHOD, not_status='DONE')) > 0:
        return
    #
    for dcd_fn in input_dcd:
        if not dcd_fn.status(): 
            continue
        run_home = dcd_fn.dirname()
        input_s = [dcd_fn]
        output_s = [run_home.fn("statpot.dat")]
        job.add_task(METHOD, input_s, output_s, use_gpu=False, n_proc=24)
    #
    job.to_json()

def run(job):
    task_s = job.get_task(METHOD, host=HOSTNAME, status='RUN') 
    if len(task_s) == 0:
        return
    #
    for index,task in task_s:
        input_dcd = task['input'][0]
        run_home = input_dcd.dirname()
        output_dat = task['output'][0]
        if output_dat.status():
            continue
        #
        run_home.chdir()
        #
        pdblist = run_home.fn("pdb_s")
        if not run_home.subdir("ens").status() or not pdblist.status():
            with pdblist.open("wt") as fout:
                cmd = ['%s/pdb_extract'%EXEC_HOME, job.top_fn.short()]
                cmd.extend(['--dir', 'ens'])
                cmd.extend(['--name', 'sample'])
                cmd.append("--structured")
                cmd.extend(['--dcd', input_dcd.short()])
                system(cmd, outfile=fout)

        cmd = [EXEC]
        cmd.extend(['-j', '%d'%task['resource'][3]])
        cmd.extend(['-l', pdblist.short()])
        cmd.append('--rwplus')
        cmd.append('--dfire')
        with output_dat.open("wt") as fout:
            system(cmd, outfile=fout, errfile='/dev/null')
    #

def submit(job):
    pass

def status(job):
    pass

def main():
    arg = argparse.ArgumentParser(prog='score')
    arg.add_argument(dest='command', choices=['prep', 'run'], help='exec type')
    arg.add_argument(dest='work_dir', help='work_dir, which has a JSON file')
    arg.add_argument('-i', '--input', dest='input_dcd', nargs='*', \
            help='input DCD file, mandatory for "prep"')  

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
    if arg.command == 'prep':
        if arg.input_dcd is None:
            sys.exit("Error: input_dcd required\n")
        arg.input_dcd = [path.Path(fn) for fn in arg.input_dcd]
        #
        prep(job, arg.input_dcd)

    elif arg.command == 'run':
        run(job)

if __name__ == '__main__':
    main()
