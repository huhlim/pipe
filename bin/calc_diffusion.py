#!/usr/bin/env python

import os
import sys
import path
import argparse

from libcommon import *

METHOD = 'calc_diffusion'
EXEC = f'{EXEC_HOME}/calc_translational_diffusion.py'

def prep(job, top_fn, dcd_fn_s, *arg, **kwarg):
    job.analysis_home = job.work_home.subdir("analysis", build=True)
    #
    input_s = [path.Path(top_fn), [path.Path(fn) for fn in dcd_fn_s]]
    output_s = [job.analysis_home.fn("translational_diffusion.pkl")]
    job.add_task(METHOD, input_s, output_s, n_proc=4, *arg, **kwarg)
    #
    job.to_json()

def run(job):
    task_s = job.get_task(METHOD, status='SUBMIT') 
    if len(task_s) == 0:
        return
    #
    job.work_home.chdir()
    #
    for index,task in task_s:
        input_s = task['input']
        output_s = task['output']
        options = task['etc']
        #
        status = True
        for output in output_s:
            if not output.status():
                status = False ; break
        if status: continue
        #
        cmd = [EXEC]
        cmd.append(output_s[0].short())
        cmd.extend(["--top", input_s[0].short()])
        cmd.extend(['--traj'] + [fn.short() for fn in input_s[1:]])
        for key,value in options.items():
            if key in ['overwrite']:
                if value: cmd.append(f"--{key}")
                continue
            else:
                cmd.append(f"--{key}")
            #
            if isinstance(value, int):
                cmd.append('%d'%value)
            elif isinstance(value, float):
                cmd.append('%f'%value)
            elif isinstance(value, str):
                cmd.append(value)
            elif isinstance(value, list):
                for v in value:
                    if isinstance(v, str):
                        cmd.append(v)
                    elif isinstance(v, path.Path):
                        cmd.append(v.short())

        system(cmd)

def submit(job):
    task_s = job.get_task(METHOD, status='SUBMIT') 
    if len(task_s) == 0:
        return
    #
    job.work_home.chdir()
    #
    for index,task in task_s:
        input_s = task['input']
        output_s = task['output']
        options = task['etc']
        #
        status = True
        for output in output_s:
            if not output.status():
                status = False ; break
        if status: continue
        #
        cmd_s = []
        cmd_s.append("cd %s\n"%(job.work_home))
        #
        cmd = [EXEC]
        cmd.append(output_s[0].short())
        cmd.extend(["--top", input_s[0].short()])
        cmd.extend(['--traj'] + [path.Path(fn).short() for fn in input_s[1]])
        #
        for key,value in options.items():
            if key in ['overwrite']:
                if value: cmd.append(f"--{key}")
                continue
            else:
                cmd.append(f"--{key}")
            #
            if isinstance(value, int):
                cmd.append('%d'%value)
            elif isinstance(value, float):
                cmd.append('%f'%value)
            elif isinstance(value, str):
                cmd.append(value)
            elif isinstance(value, list):
                for v in value:
                    if isinstance(v, str):
                        cmd.append(v)
                    elif isinstance(v, path.Path):
                        cmd.append(v.short())

        cmd_s.append(" ".join(cmd) + '\n')
        #
        job.write_submit_script(METHOD, index, cmd_s)

