#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

METHOD = 'average_meta'
EXEC = '%s/average_meta.py'%EXEC_HOME

PARAM = {}
PARAM['score'] = ("RWplus", 25.0)
PARAM['casp12'] = ("RWplus", 0.5, 225., 45.)
PARAM['cluster'] = (2.0, 20, 5) # rmsd_cutoff, subsample, max_run_md

def prep(job, output_prefix, input_prod, input_json, rule='score'):
#    if len(job.get_task(METHOD, not_status='DONE')) > 0:
#        return
    #
    prod_s = []
    for _,prod in job.get_task("prod_meta"):
        prod_index = int(prod['output'][0].dirname().dirname().name())
        if prod_index in input_prod:
            prod_s.append(prod)
    #
    score_s = job.get_task("score")
    #
    job.average_home = job.work_home.subdir("average", build=True)
    job.average_home.chdir()
    #
    with input_json.open() as fp:
        options = json.load(fp)
    options['rule'] = (rule, PARAM[rule])
    options['input'] = {}
    options['input']['topIndex'] = []
    options['input']['top_fn'] = []
    options['input']['dcd_fn'] = []
    if rule in ['score', 'casp12']:
        options['input']['score_fn'] = []
    if rule in ['casp12']:
        options['input']['qual_fn'] = []
    #
    for prod in prod_s:
        if prod['resource'][0] != 'DONE':
            return
        if not prod['output'][0].status():
            return

        if rule in ['score', 'casp12']:
            score = None
            for _,s in score_s:
                if s['input'][0] == prod['output'][0]:
                    score = s
                    break
            if score is None:
                return
            if score['resource'][0] != 'DONE':
                return
            if not score['output'][0].status():
                return
        #
        options['input']['topIndex'].append(prod['input'][1].path())
        options['input']['top_fn'].append(prod['input'][2].path())
        options['input']['dcd_fn'].append(prod['output'][0].path())
        if rule in ['score', 'casp12']:
            options['input']['score_fn'].append(score['output'][0].path())
        if rule in ['casp12']:
            options['input']['qual_fn'].append(score['output'][1].path())
    #
    job_input_json = job.average_home.fn("%s.input.json"%output_prefix)
    with job_input_json.open("wt") as fout:
        fout.write(json.dumps(options, indent=2, default=JSONserialize))
    #
    input_s = [output_prefix, job_input_json]
    output_s = [job.average_home.fn("%s.pdb_s"%output_prefix)]
    job.add_task(METHOD, input_s, output_s, use_gpu=True, n_proc=1)
    #
    job.to_json()

def run(job):
    task_s = job.get_task(METHOD, host=HOSTNAME, status='RUN') 
    if len(task_s) == 0:
        return
    gpu_id = os.environ['CUDA_VISIBLE_DEVICES']
    #
    for index,task in task_s:
        if task['resource'][1].split("/")[1] != gpu_id: continue
        input_s = task['input']
        output_prefix = input_s[0]
        input_json = input_s[1]
        output_pdb = task['output'][0]
        if output_pdb.status():
            continue
        #
        #with input_json.open() as fp:
        #    options = json.load(fp)
        #options['ssbond'] = []
        #for line in job.ssbond:
        #    chain_1 = line[15]
        #    chain_2 = line[29]
        #    if chain_1 == ' ' and chain_2 == ' ':
        #        line = '%sA%sA%s'%(line[:15], line[16:29], line[30:])
        #    options['ssbond'].append(line)
        #
        job.average_home.chdir()
        #
        #input_json = job.average_home.fn("%s.json"%output_prefix)
        #with input_json.open("wt") as fout:
        #    fout.write(json.dumps(options, indent=2, default=JSONserialize))
        #
        cmd = [EXEC, output_prefix]
        cmd.extend(['--input', input_json.short()])
        #
        system(cmd, verbose=job.verbose)

def submit(job):
    task_s = job.get_task(METHOD, status='SUBMIT')
    if len(task_s) == 0:
        return
    #
    for index,task in task_s:
        if task['resource'][1].split("/")[1] != gpu_id: continue
        input_s = task['input']
        output_prefix = input_s[0]
        input_json = input_s[1]
        output_pdb = task['output'][0]
        if output_pdb.status():
            continue
        #
        #with input_json.open() as fp:
        #    options = json.load(fp)
        #options['ssbond'] = []
        #for line in job.ssbond:
        #    chain_1 = line[15]
        #    chain_2 = line[29]
        #    if chain_1 == ' ' and chain_2 == ' ':
        #        line = '%sA%sA%s'%(line[:15], line[16:29], line[30:])
        #    options['ssbond'].append(line)
        #options['rule'] = rule
        #
        job.average_home.chdir()
        #
        #input_json = job.average_home.fn("%s.json"%output_prefix)
        #with input_json.open("wt") as fout:
        #    fout.write(json.dumps(options, indent=2, default=JSONserialize))
        #
        cmd_s = []
        cmd_s.append("cd %s\n"%average_home)
        cmd = [EXEC, output_prefix]
        cmd.extend(['--input', input_json.short()])
        cmd_s.append(" ".join(cmd) + '\n')
        #
        job.write_submit_script(METHOD, index, cmd_s)

def status(job):
    pass

def main():
    arg = argparse.ArgumentParser(prog='average')
    arg.add_argument(dest='command', choices=['prep', 'run'], help='exec type')
    arg.add_argument(dest='work_dir', help='work_dir, which has a JSON file')
    arg.add_argument('-o', '--output', dest='output_prefix',\
            help='output prefix, mandatory for "prep"')
    arg.add_argument('-i', '--input', dest='input_prod', type=int, nargs='*', \
            help='input PRODs, mandatory for "prep"')  

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
        if arg.input_prod is None:
            sys.exit("Error: input_prod required\n")
        #
        prep(job, arg.output_prefix, arg.input_prod)

    elif arg.command == 'run':
        run(job)

if __name__ == '__main__':
    main()
