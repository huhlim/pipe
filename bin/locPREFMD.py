#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

METHOD = 'locPREFMD'
EXEC = '%s/locprefmd.sh'%EXEC_HOME

def prep(job, input_pdb):
    if len(job.get_task(METHOD, not_status='DONE')) > 0:
        return
    #
    for fn in input_pdb:
        out = fn.dirname().fn("%s.prefmd.pdb"%(fn.name()))
        #
        if CHARMM_MPI:
            job.add_task(METHOD, [fn], [out], use_gpu=False, n_proc=16)
        else:
            job.add_task(METHOD, [fn], [out], use_gpu=False, n_proc=4)
    #
    job.to_json()

def run(job):
    task_s = job.get_task(METHOD, host=HOSTNAME, status='RUN') 
    if len(task_s) == 0:
        return
    os.environ['CHARMMEXEC'] = CHARMMEXEC_MPI
    #
    for index,task in task_s:
        input_pdb = task['input'][0]
        run_home = input_pdb.dirname()
        output_pdb = task['output'][0]
        if output_pdb.status():
            continue
        #
        run_home.chdir()
        with output_pdb.open("wt") as fout:
            cmd = [EXEC, input_pdb.short()]
            system(cmd, outfile=fout, verbose=job.verbose)

def submit(job):
    task_s = job.get_task(METHOD, status='WAIT') 
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
        cmd = []
        cmd.append("export CHARMMEXEC=%s\n\n"%CHARMMEXEC_MPI)
        cmd.append("cd %s\n"%run_home)
        cmd.append("%s %s > %s\n"%(EXEC, input_pdb.short(), output_pdb.short()))
        #
        job.write_submit_script(METHOD, index, cmd, submit=True)

def status(job):
    task_s = job.get_task(METHOD)
    if len(task_s) == 0:
        return
    #
    for index,task in task_s:
        output_pdb = task['output'][0]
        if output_pdb.status():
            job.update_task_status(METHOD, index, "DONE")
        elif output_pdb.exists():
            job.update_task_status(METHOD, index, "RUN")

def main():
    arg = argparse.ArgumentParser(prog='locPREFMD')
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
