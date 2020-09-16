#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

METHOD = 'equil'
EXEC_SOLUBLE = '%s/equil.py'%EXEC_HOME
EXEC_MEMBRANE = '%s/equil_membrane.py'%EXEC_HOME

def prep(job, equil_index, input_pdb, input_json):
#    if len(job.get_task(METHOD, not_status='DONE')) > 0:
#        return
    #
    OUTs = ['%s.psf'%job.title, '%s.orient.pdb'%job.title, '%s.equil.restart.pkl'%job.title, '%s.equil.pdb'%job.title]
    #
    job.equil_home = job.work_home.subdir("equil", build=True)
    job.equil_home.chdir()
    #
    for fn in input_pdb:
        run_home = job.equil_home.subdir("%d"%equil_index)
        #
        input_s = [run_home, fn, input_json]
        output_s = [run_home.fn(X) for X in OUTs]
        status = True
        for output in output_s:
            if not output.status():
                status = False ; break
        #
        if status: 
            job.add_task(METHOD, input_s, output_s, use_gpu=True, n_proc=12, status='DONE')
        else:
            job.add_task(METHOD, input_s, output_s, use_gpu=True, n_proc=12)
        equil_index += 1
    #
    job.to_json()

def prep_membrane(job, equil_index, input_pdb, input_psf, input_crd, input_json):
    if len(job.get_task(METHOD, not_status='DONE')) > 0:
        return
    #
    OUTs = ['%s.psf'%job.title, '%s.orient.pdb'%job.title, '%s.equil.restart.pkl'%job.title, '%s.equil.pdb'%job.title]
    #
    job.equil_home = job.work_home.subdir("equil", build=True)
    job.equil_home.chdir()
    #
    for pdb_fn, psf_fn, crd_fn in zip(input_pdb, input_psf, input_crd):
        run_home = job.equil_home.subdir("%d"%equil_index)
        #
        input_s = [run_home, pdb_fn, psf_fn, crd_fn, input_json]
        output_s = [run_home.fn(X) for X in OUTs]
        status = True
        for output in output_s:
            if not output.status():
                status = False ; break
        #
        if status: 
            job.add_task(METHOD, input_s, output_s, use_gpu=True, n_proc=12, status='DONE')
        else:
            job.add_task(METHOD, input_s, output_s, use_gpu=True, n_proc=12)
        equil_index += 1
    #
    job.to_json()

def run(job):
    task_s = job.get_task(METHOD, host=HOSTNAME, status='RUN') 
    if len(task_s) == 0:
        return
    os.environ['CHARMMEXEC'] = CHARMMEXEC
    gpu_id = os.environ['CUDA_VISIBLE_DEVICES']
    #
    for index,task in task_s:
        if task['resource'][1].split("/")[1] != gpu_id: continue
        #
        run_home = task['input'][0]
        input_pdb  = task['input'][1]
        input_json = task['input'][-1]
        #
        if job.has("is_membrane_protein"):
            is_membrane = True
            EXEC = EXEC_MEMBRANE
            input_psf = task['input'][2]
            input_crd = task['input'][3]
        else:
            is_membrane = False
            EXEC = EXEC_SOLUBLE
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
        options['ssbond'] = []
        for line in job.ssbond:
            chain_1 = line[15]
            chain_2 = line[29]
            if chain_1 == ' ' and chain_2 == ' ':
                line = '%sA%sA%s'%(line[:15], line[16:29], line[30:])
            options['ssbond'].append(line)
        #
        options['input_pdb'] = input_pdb
        options['input_json'] = input_json
        if job.has("has_ligand"):
            options['ligand_json'] = job.ligand_json
        #
        run_home.build()
        run_home.chdir()
        #
        equil_json = run_home.fn("input.json")
        if not equil_json.status():
            with equil_json.open("wt") as fout:
                fout.write(json.dumps(options, indent=2, default=JSONserialize))
        #
        if not is_membrane:
            cmd = [EXEC, job.title, input_pdb.short()]
        else:
            cmd = [EXEC, job.title, input_pdb.short(), input_psf.short(), input_crd.short()]
        cmd.extend(['--input', equil_json.short()])
        if job.verbose:  cmd.append('--verbose')
        if job.keep_tmp: cmd.append('--keep')
        #
        system(cmd, verbose=job.verbose)

def submit(job):
    task_s = job.get_task(METHOD, status='SUBMIT') 
    if len(task_s) == 0:
        return
    #
    for index,task in task_s:
        run_home = task['input'][0]
        input_pdb  = task['input'][1]
        input_json = task['input'][-1]
        #
        if job.has("is_membrane_protein"):
            is_membrane = True
            EXEC = EXEC_MEMBRANE
            input_psf = task['input'][2]
            input_crd = task['input'][3]
        else:
            is_membrane = False
            EXEC = EXEC_SOLUBLE
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
        options['ssbond'] = []
        for line in job.ssbond:
            chain_1 = line[15]
            chain_2 = line[29]
            if chain_1 == ' ' and chain_2 == ' ':
                line = '%sA%sA%s'%(line[:15], line[16:29], line[30:])
            options['ssbond'].append(line)
        #
        options['input_pdb'] = input_pdb
        options['input_json'] = input_json
        if job.has("has_ligand"):
            options['ligand_json'] = job.ligand_json
        #
        run_home.build()
        run_home.chdir()
        #
        equil_json = run_home.fn("input.json")
        if not equil_json.status():
            with equil_json.open("wt") as fout:
                fout.write(json.dumps(options, indent=2, default=JSONserialize))
        #
        cmd_s = []
        cmd_s.append("cd %s\n"%run_home)
        if not is_membrane:
            cmd = [EXEC, job.title, input_pdb.short()]
        else:
            cmd = [EXEC, job.title, input_pdb.short(), input_psf.short(), input_crd.short()]
        cmd.extend(['--input', equil_json.short()])
        if job.verbose:  cmd.append('--verbose')
        if job.keep_tmp: cmd.append('--keep')
        cmd_s.append(" ".join(cmd) + '\n')
        #
        job.write_submit_script(METHOD, index, cmd_s)

def status(job):
    task_s = job.get_task(METHOD)
    if len(task_s) == 0:
        return
    #
    for index,task in task_s:
        run_home = task['input'][0]
        output_s = task['output']
        status = True
        for output in output_s:
            if not output.status():
                status = False ; break
        #
        if status:
            job.update_task_status(METHOD, index, "DONE")
        elif run_home.exists():
            job.update_task_status(METHOD, index, "RUN")

def main():
    arg = argparse.ArgumentParser(prog='equil')
    arg.add_argument(dest='command', choices=['prep', 'run'], help='exec type')
    arg.add_argument(dest='work_dir', help='work_dir, which has a JSON file')
    arg.add_argument('--index', dest='equil_index', help='equil_index', type=int)
    arg.add_argument('-i', '--input', dest='input_pdb', nargs='*', \
            help='input PDB file, mandatory for "prep"')  
    arg.add_argument('-j', '--json', dest='input_json', \
            help='input JSON file, mandatory for "prep"')  

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
        prep(job, arg.equil_index, arg.input_pdb, path.Path(arg.input_json))

    elif arg.command == 'run':
        run(job)

if __name__ == '__main__':
    main()
