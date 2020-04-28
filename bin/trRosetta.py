#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

METHOD = 'trRosetta'
EXEC = '%s/trRosetta/predict.py'%EXEC_HOME

def prep(job, input_fa):
    if len(job.get_task(METHOD, not_status='DONE')) > 0:
        return
    #
    job.trRosetta_home = job.work_home.subdir("trRosetta", build=True)
    out = job.trRosetta_home.fn("model_s")
    #
    job.add_task(METHOD, [job.title, input_fa, job.trRosetta_home], [out], use_gpu=False, n_proc=48)
    #
    job.to_json()

def run(job):
    task_s = job.get_task(METHOD, host=HOSTNAME, status='RUN') 
    if len(task_s) == 0:
        return
    #
    for index,task in task_s:
        title = task['input'][0]
        input_fa = task['input'][1]
        run_home = task['input'][2]
        output_list = task['output'][0]
        if output_list.status():
            continue
        #
        run_home.chdir()
        cmd = [EXEC, title, input_fa.short()]
        system(cmd, verbose=job.verbose)

def submit(job):
    task_s = job.get_task(METHOD, status='SUBMIT') 
    if len(task_s) == 0:
        return
    #
    for index,task in task_s:
        title = task['input'][0]
        input_fa = task['input'][1]
        run_home = task['input'][2]
        output_list = task['output'][0]
        if output_list.status():
            continue
        #
        run_home.chdir()
        #
        cmd = []
        cmd.append("cd %s\n"%run_home)
        cmd.append(" ".join([EXEC, title, input_fa.short()]) + '\n')
        #
        job.write_submit_script(METHOD, index, cmd)

def main():
    arg = argparse.ArgumentParser(prog='trRosetta')
    arg.add_argument(dest='command', choices=['prep', 'run'], help='exec type')
    arg.add_argument(dest='work_dir', help='work_dir, which has a JSON file')
    arg.add_argument('-i', '--input', dest='input_fa', \
            help='input FASTA file, mandatory for "prep"')  

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
        if arg.input_fa is None:
            sys.exit("Error: input_fa required\n")
        arg.input_fa = path.Path(arg.input_fa)
        #
        prep(job, arg.input_fa)

    elif arg.command == 'run':
        run(job)

if __name__ == '__main__':
    main()
