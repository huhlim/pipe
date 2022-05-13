#!/usr/bin/env python

import os
import sys
import path
import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt

from libcommon import *

sys.path.append('%s/bin/exec'%os.getenv("PIPE_HOME"))
from libanalysis import get_molecules

METHOD = 'diffusion'
EXEC = f'{EXEC_HOME}/calc_translational_diffusion.py'

def prep(job, top_fn, dcd_fn_s, *arg, **kwarg):
    job.analysis_home = job.work_home.subdir("analysis", build=True)
    #
    for dcd_fn in dcd_fn_s:
        run_name = '/'.join(dcd_fn.dirname().split("/")[-2:])
        run_home = job.analysis_home.subdir(run_name, build=True)
        #
        input_s = [path.Path(top_fn), path.Path(dcd_fn)]
        output_s = [run_home.fn("translational_diffusion.pkl")]
        job.add_task(METHOD, input_s, output_s, n_proc=4, *arg, **kwarg)
    #
    job.to_json()

def read_pkl(pkl_fn_s, group_s):
    _msd_s = {group: {} for group in group_s}
    _diffusion_s = {group: {} for group in group_s}
    #
    for pkl_fn in pkl_fn_s:
        with pkl_fn.open("rb") as fp:
            X = pickle.load(fp)
        msd = X['msd']
        diffusion = X['diffusion']
        #
        for group, name_s in group_s.items():
            for name in name_s:
                if name not in msd:
                    raise KeyError(name)
                for k in range(len(msd[name])): # number of traj 
                    key = (pkl_fn.short(), k)
                    if key not in _msd_s[group]:
                        _msd_s[group][key] = []
                        _diffusion_s[group][key] = []
                    #
                    t_lag = msd[name][k][0]
                    #
                    _msd_s[group][key].append(msd[name][k][1])
                    _diffusion_s[group][key].append(diffusion[name][k])
    #
    msd_s = {group: [] for group in group_s}
    diffusion_s = {group: [] for group in group_s}
    for group in group_s:
        for runner in _msd_s[group]:
            msd = np.concatenate(_msd_s[group][runner], axis=0).mean(axis=0)
            diffusion = np.concatenate(_diffusion_s[group][runner], axis=0).mean()
            msd_s[group].append(msd)
            diffusion_s[group].append(diffusion)
    return t_lag, msd_s, diffusion_s

def plot(out_home, t_lag, msd_s):
    for group,msd in msd_s.items():
        png_fn = out_home.fn(f"translational_diffusion.{group}.png")
        #
        m = np.mean(msd, 0)
        s = np.std(msd, 0) / np.sqrt(len(msd))
        #
        fig, ax = plt.subplots(figsize=(4.8, 4.8))
        #
        ax.plot(t_lag, m, 'r-', linewidth=1.5)
        ax.fill_between(t_lag, m-s, m+s, color='red', alpha=0.2)
        #
        ax.set_xlim((0, 100))
        ax.set_xticks(np.arange(0, 101, 20))
        ax.set_ylim(bottom=0)
        #
        fig.tight_layout()
        #
        plt.savefig(png_fn.short())
        plt.close("all")

def summarize(job, pkl_fn_s, group_s):
    out_fn = job.analysis_home.fn("translational_diffusion.summary.pkl")
    status = out_fn.status()
    for group in group_s:
        png_fn = job.analysis_home.fn(f"translational_diffusion.{group}.png")
        if not png_fn.status():
            status = False ; break
    if status: return
    #
    t_lag, msd_s, diffusion_s = read_pkl(pkl_fn_s, group_s)
    plot(job.analysis_home, t_lag, msd_s)
    #
    with out_fn.open("wb") as fout:
        pickle.dump({"msd": (t_lag, msd_s), "diffusion": diffusion_s}, fout)

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
        cmd.extend(['--traj', input_s[1].short()])
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
        cmd.extend(['--traj', input_s[1].short()])
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

