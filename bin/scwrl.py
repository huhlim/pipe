#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

METHOD = 'scwrl'
EXEC = 'scwrl4'

def prep(job, input_pdb):
    if len(job.get_task(METHOD, not_status='DONE')) > 0:
        return
    #
    for pdb_fn in input_pdb:
        out_fn = pdb_fn.dirname().fn("%s.scwrl.pdb"%(pdb_fn.name()))
        if out_fn.status(): 
            continue
        input_s = [pdb_fn]
        output_s = [out_fn]
        job.add_task(METHOD, input_s, output_s, use_gpu=False, n_proc=1)
    #
    job.to_json()

def run(job):
    task_s = job.get_task(METHOD, host=HOSTNAME, status='RUN') 
    if len(task_s) == 0:
        return
    #
    for index,task in task_s:
        input_pdb = task['input'][0]
        run_home = input_pdb.dirname()
        output_pdb = task['output'][0]
        if output_pdb.status():
            continue
        #
        run_home.chdir()
        #
        cmd = [EXEC]
        cmd.extend(['-i', input_pdb.short()])
        cmd.extend(['-o', output_pdb.short()])
        #
        system(cmd, stdout=True)

def submit(job):
    task_s = job.get_task(METHOD, host=HOSTNAME, status='RUN') 
    if len(task_s) == 0:
        return
    #
    que = []
    for index,task in task_s:
        input_pdb = task['input'][0]
        run_home = input_pdb.dirname()
        output_pdb = task['output'][0]
        if output_pdb.status():
            continue
        #
        run_home.chdir()
        #
        que.append("cd %s\n"%run_home)
        #
        cmd = [EXEC]
        cmd.extend(['-i', input_pdb.short()])
        cmd.extend(['-o', output_pdb.short()])
        cmd.append("&> /dev/null\n")
        que.append(" ".join(cmd))
    #
    job.write_submit_script(METHOD, 0, que)

def status(job):
    pass

def main():
    arg = argparse.ArgumentParser(prog='scwrl')
    arg.add_argument(dest='command', choices=['prep', 'run'], help='exec type')
    arg.add_argument(dest='work_dir', help='work_dir, which has a JSON file')
    arg.add_argument('-i', '--input', dest='input_pdb', nargs='*', \
            help='input PDB file, mandatory for "prep"')  

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
        if arg.input_pdb is None:
            sys.exit("Error: input_pdb required\n")
        arg.input_pdb = [path.Path(fn) for fn in arg.input_pdb]
        #
        prep(job, arg.input_pdb)

    elif arg.command == 'run':
        run(job)

if __name__ == '__main__':
    main()
