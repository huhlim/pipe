#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

METHOD = 'prod'
EXEC = '%s/prod.py'%EXEC_HOME

def prep(job, prod_index, input_equil, input_json, n_replica):
    if len(job.get_task(METHOD, not_status='DONE')) > 0:
        return
    #
    job.prod_home = job.work_home.subdir("prod", build=True)
    job.prod_home.chdir()
    #
    iter_home = job.prod_home.subdir("%d"%prod_index, build=True)
    options = {}
    options['input_equil'] = input_equil
    options['input_json'] = input_json.path()
    options['n_replica'] = n_replica
    with iter_home.fn("input.json").open("wt") as fout:
        fout.write(json.dumps(options, indent=2))
    #
    for i in range(n_replica):
        run_home = iter_home.subdir("%d"%i)
        #
        input_s = [run_home, input_equil, input_json]
        output_s = [run_home.fn("solute.dcd")]
        job.add_task(METHOD, input_s, output_s, use_gpu=True, n_proc=1)
    #
    job.to_json()

def run(job):
    task_s = job.get_task(METHOD, host=HOSTNAME, status='RUN') 
    if len(task_s) == 0:
        return
    #
    for index,task in task_s:
        run_home = task['input'][0]
        input_equil_index  = task['input'][1]
        input_json = task['input'][2]
        #
        equil_home = job.get_task("equil")[input_equil_index][1]['input'][0]
        #
        output_s = task['output']
        status = True
        for output in output_s:
            if not output.status():
                status = False ; break
        if status: continue
        #
        with input_json.open() as fp:
            options = json.load(fp)
        #
        run_home.build()
        run_home.chdir()
        run_name = 'r%s'%(run_home.split("/")[-1])
        #
        if 'restraint' in options:
            options['restraint']['reference'] = equil_home.fn("%s.orient.pdb"%job.title).short()
        options['input'] = {}
        options['input']['psf'] = equil_home.fn("%s.psf"%job.title).short()
        options['input']['pdb'] = equil_home.fn("%s.equil.pdb"%job.title).short()
        options['input']['n_atom'] = job.n_atom
        if options['restart']:
            options['input']['restart'] = equil_home.fn("%s.equil.restart.pkl"%job.title).short()
        #
        run_json = run_home.fn("input.json")
        if not run_json.status():
            with run_json.open("wt") as fout:
                fout.write(json.dumps(options, indent=2))
        #
        cmd = [EXEC, run_name]
        cmd.extend(["--input", run_json.short()])
        if job.verbose:  cmd.append('--verbose')
        if job.keep_tmp: cmd.append('--keep')
        #
        system(cmd, verbose=job.verbose)

def submit(job):
    pass

def status(job):
    pass

def main():
    arg = argparse.ArgumentParser(prog='prod')
    arg.add_argument(dest='command', choices=['prep', 'run'], help='exec type')
    arg.add_argument(dest='work_dir', help='work_dir, which has a JSON file')
    arg.add_argument('--index', dest='prod_index', help='prod_index', type=int)
    arg.add_argument('-i', '--input', dest='input_equil', nargs='*', \
            help='input equil index, mandatory for "prep"')  
    arg.add_argument('-j', '--json', dest='input_json', \
            help='input JSON file, mandatory for "prep"')  
    arg.add_argument('-n', dest='n_replica', type=int, default=5, \
            help='number of replicas')

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
        if arg.input_json is None:
            sys.exit("Error: input_json required\n")
        arg.input_pdb = [path.Path(fn) for fn in arg.input_pdb]
        #
        prep(job, arg.prod_index, arg.input_equil, path.Path(arg.input_json), arg.n_replica)

    elif arg.command == 'run':
        run(job)

if __name__ == '__main__':
    main()
