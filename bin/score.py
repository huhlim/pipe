#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

METHOD = 'score'
EXEC1 = '%s/calc_statpot'%EXEC_HOME
EXEC2 = '%s/calc_tmscore'%EXEC_HOME

def prep(job, input_dcd):
    if len(job.get_task(METHOD, not_status='DONE')) > 0:
        return
    #
    for dcd_fn in input_dcd:
        if not dcd_fn.status(): 
            continue
        run_home = dcd_fn.dirname()
        input_s = [dcd_fn]
        output_s = [run_home.fn("statpot.dat"), run_home.fn("qual_init.dat")]
        job.add_task(METHOD, input_s, output_s, use_gpu=False, n_proc=16)
    #
    job.to_json()

def prep_meta(job, prod_task_s):
    if len(job.get_task(METHOD, not_status='DONE')) > 0:
        return
    #
    for _,task in prod_task_s:
        if task['resource'][0] != 'DONE':
            continue
        #
        dcd_fn = task['output'][0]
        run_home = task['input'][0]
        index_fn = task['input'][1]
        top_fn = task['input'][2]
        #
        input_s = [dcd_fn]
        output_s = [run_home.fn("statpot.dat"), run_home.fn("qual_init.dat")]
        etc_s = [index_fn, top_fn]
        #
        job.add_task(METHOD, input_s, output_s, \
                index_fn=index_fn, top_fn=top_fn, \
                use_gpu=False, n_proc=16)
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
        output_score = task['output'][0]
        output_qual = task['output'][1]
        if output_score.status() and output_qual.status():
            continue
        #
        index_fn = task['etc'].get("index_fn", None)
        top_fn = task['etc'].get("top_fn", job.top_fn)
        #if len(task['etc']) > 0:    # meta
        #    index_fn = task['etc'][0]
        #    top_fn = task['etc'][1]
        #else:
        #    index_fn = None
        #    top_fn = job.top_fn
        #
        run_home.chdir()
        #
        pdblist = run_home.fn("pdb_s")
        if not run_home.subdir("ens").status() or not pdblist.status():
            with pdblist.open("wt") as fout:
                cmd = ['%s/pdb_extract'%EXEC_HOME, top_fn.short()]
                if index_fn is not None:
                    cmd.extend(['--topIndex', index_fn.short()])
                cmd.extend(['--dir', 'ens'])
                cmd.extend(['--name', 'sample'])
                cmd.append("--structured")
                cmd.extend(['--dcd', input_dcd.short()])
                system(cmd, outfile=fout)

        if not output_score.status():
            cmd = [EXEC1]
            cmd.extend(['-j', '%d'%task['resource'][3]])
            cmd.extend(['-l', pdblist.short()])
            cmd.append('--rwplus')
            cmd.append('--dfire')
            with output_score.open("wt") as fout:
                system(cmd, outfile=fout, errfile='/dev/null')

        if not output_qual.status():
            cmd = [EXEC2]
            cmd.extend(['-j', '%d'%task['resource'][3]])
            cmd.extend(['-r', job.init_pdb[0].short()])
            cmd.extend(['-l', pdblist.short()])
            with output_qual.open("wt") as fout:
                system(cmd, outfile=fout, errfile='/dev/null')
    #

def submit(job):
    task_s = job.get_task(METHOD, status='SUBMIT') 
    if len(task_s) == 0:
        return
    #
    for index,task in task_s:
        input_dcd = task['input'][0]
        run_home = input_dcd.dirname()
        output_score = task['output'][0]
        output_qual = task['output'][1]
        if output_score.status() and output_qual.status():
            continue
        #
        cmd_s = []
        run_home.chdir()
        cmd_s.append("cd %s\n"%run_home)
        #
        pdblist = run_home.fn("pdb_s")
        if not run_home.subdir("ens").status() or not pdblist.status():
            cmd = ['%s/pdb_extract'%EXEC_HOME, job.top_fn.short()]
            cmd.extend(['--dir', 'ens'])
            cmd.extend(['--name', 'sample'])
            cmd.append("--structured")
            cmd.extend(['--dcd', input_dcd.short()])
            cmd.append("> %s 2> /dev/null"%pdblist.short())
            cmd_s.append(' '.join(cmd)+'\n')

        if not output_score.status():
            cmd = [EXEC1]
            cmd.extend(['-j', '%d'%task['resource'][3]])
            cmd.extend(['-l', pdblist.short()])
            cmd.append('--rwplus')
            cmd.append('--dfire')
            cmd.append("> %s 2> /dev/null"%output_score.short())
            cmd_s.append(' '.join(cmd)+'\n')

        if not output_qual.status():
            cmd = [EXEC2]
            cmd.extend(['-j', '%d'%task['resource'][3]])
            cmd.extend(['-r', job.init_pdb[0].short()])
            cmd.extend(['-l', pdblist.short()])
            cmd.append("> %s 2> /dev/null"%output_qual.short())
            cmd_s.append(' '.join(cmd)+'\n')
        #
        job.write_submit_script(METHOD, index, cmd_s)

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
