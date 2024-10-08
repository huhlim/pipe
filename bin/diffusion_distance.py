#!/usr/bin/env python

import os
import sys
import path
import argparse
import pickle
import numpy as np
import itertools
import matplotlib.pyplot as plt

from libcommon import *

METHOD = "diffusion_distance"
EXEC = f"{EXEC_HOME}/calc_distance_diffusion.py"


def prep(job, top_fn, dcd_fn_s, *arg, **kwarg):
    job.analysis_home = job.work_home.subdir("analysis", build=True)
    #
    for dcd_fn in dcd_fn_s:
        run_name = "/".join(dcd_fn.dirname().split("/")[-2:])
        run_home = job.analysis_home.subdir(run_name, build=True)
        #
        input_s = [path.Path(top_fn), path.Path(dcd_fn)]
        output_s = [run_home.fn("distance_diffusion.pkl")]
        job.add_task(METHOD, input_s, output_s, n_proc=8, *arg, **kwarg)
    #
    job.to_json()


def read_pkl(pkl_fn_s, group_s):
    _msd_s = {}
    _diffusion_s = {}
    #
    for pkl_fn in pkl_fn_s:
        with pkl_fn.open("rb") as fp:
            X = pickle.load(fp)
        msd = X["msd"]
        diffusion = X["diffusion"]
        #
        for (group_i, group_j) in itertools.combinations_with_replacement(group_s, 2):
            group_pair = (group_i, group_j)
            name_i = group_s[group_i]
            name_j = group_s[group_j]
            #
            pair_s = itertools.product(name_i, name_j)
            processed = []
            for pair in pair_s:
                if pair not in msd:
                    pair = (pair[1], pair[0])
                if pair not in msd:
                    continue
                if pair in processed:
                    continue
                processed.append(pair)
                #
                if group_pair not in _msd_s:
                    _msd_s[group_pair] = {}
                    _diffusion_s[group_pair] = {}
                #
                for k in range(len(msd[pair])):  # number of traj
                    key = (pkl_fn.short(), k)
                    if key not in _msd_s[group_pair]:
                        _msd_s[group_pair][key] = []
                        _diffusion_s[group_pair][key] = []
                    #
                    t_lag = msd[pair][k][0]
                    #
                    _msd_s[group_pair][key].append(msd[pair][k][1])
                    _diffusion_s[group_pair][key].append(diffusion[pair][k])
    #
    msd_s = {}
    diffusion_s = {}
    for group_pair in _msd_s:
        msd_s[group_pair] = []
        diffusion_s[group_pair] = []
        #
        for runner in _msd_s[group_pair]:
            msd = np.array(_msd_s[group_pair][runner]).mean(axis=(0, 2))
            diffusion = np.array(_diffusion_s[group_pair][runner]).mean(axis=(0, 2))
            msd_s[group_pair].append(msd)
            diffusion_s[group_pair].append(diffusion)
    return t_lag, msd_s, diffusion_s


def plot(out_home, prefix, t_lag, msd_s):
    for group, msd in msd_s.items():
        png_fn = out_home.fn(f"distance_diffusion.{prefix}.{group[0]}-{group[1]}.png")
        #
        m = np.mean(msd, 0)
        s = np.std(msd, 0) / np.sqrt(len(msd))
        #
        fig, ax = plt.subplots(figsize=(4.8, 4.8))
        #
        ax.plot(t_lag, m, "r-", linewidth=1.5)
        ax.fill_between(t_lag, m - s, m + s, color="red", alpha=0.2)
        #
        ax.set_xlim((0, t_lag[-1]))
        ax.set_ylim(bottom=0)
        #
        fig.tight_layout()
        #
        plt.savefig(png_fn.short())
        plt.close("all")


def summarize(job, prefix, pkl_fn_s, group_s):
    out_fn = job.analysis_home.fn(f"distance_diffusion.{prefix}.summary.pkl")
    if out_fn.status():
        return
    #
    t_lag, msd_s, diffusion_s = read_pkl(pkl_fn_s, group_s)
    plot(job.analysis_home, prefix, t_lag, msd_s)
    #
    with out_fn.open("wb") as fout:
        pickle.dump({"msd": (t_lag, msd_s), "diffusion": diffusion_s}, fout)


def run(job):
    task_s = job.get_task(METHOD, status="SUBMIT")
    if len(task_s) == 0:
        return
    #
    job.work_home.chdir()
    #
    for index, task in task_s:
        input_s = task["input"]
        output_s = task["output"]
        options = task["etc"]
        #
        status = True
        for output in output_s:
            if not output.status():
                status = False
                break
        if status:
            continue
        #
        cmd = [EXEC]
        cmd.append(output_s[0].short())
        cmd.extend(["--top", input_s[0].short()])
        cmd.extend(["--traj", input_s[1].short()])
        for key, value in options.items():
            if key in ["overwrite"]:
                if value:
                    cmd.append(f"--{key}")
                continue
            else:
                cmd.append(f"--{key}")
            #
            if isinstance(value, int):
                cmd.append("%d" % value)
            elif isinstance(value, float):
                cmd.append("%f" % value)
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
    task_s = job.get_task(METHOD, status="SUBMIT")
    if len(task_s) == 0:
        return
    #
    job.work_home.chdir()
    #
    for index, task in task_s:
        input_s = task["input"]
        output_s = task["output"]
        options = task["etc"]
        #
        status = True
        for output in output_s:
            if not output.status():
                status = False
                break
        if status:
            continue
        #
        cmd_s = []
        cmd_s.append("cd %s\n" % (job.work_home))
        #
        cmd = [EXEC]
        cmd.append(output_s[0].short())
        cmd.extend(["--top", input_s[0].short()])
        cmd.extend(["--traj", input_s[1].short()])
        #
        for key, value in options.items():
            if key in ["overwrite"]:
                if value:
                    cmd.append(f"--{key}")
                continue
            else:
                cmd.append(f"--{key}")
            #
            if isinstance(value, int):
                cmd.append("%d" % value)
            elif isinstance(value, float):
                cmd.append("%f" % value)
            elif isinstance(value, str):
                cmd.append(value)
            elif isinstance(value, list):
                for v in value:
                    if isinstance(v, str):
                        cmd.append(v)
                    elif isinstance(v, path.Path):
                        cmd.append(v.short())

        cmd_s.append(" ".join(cmd) + "\n")
        #
        job.write_submit_script(METHOD, index, cmd_s)
