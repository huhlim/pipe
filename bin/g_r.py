#!/usr/bin/env python

import os
import sys
import path
import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt

from libcommon import *

METHOD = 'g_r'
EXEC = f'{EXEC_HOME}/calc_g_r.py'

def prep(job, top_fn, dcd_fn_s, *arg, **kwarg):
    job.analysis_home = job.work_home.subdir("analysis", build=True)
    n_proc=kwarg.get("n_proc", 4)
    #
    for dcd_fn in dcd_fn_s:
        run_name = '/'.join(dcd_fn.dirname().split("/")[-2:])
        run_home = job.analysis_home.subdir(run_name, build=True)
        #
        input_s = [path.Path(top_fn), path.Path(dcd_fn)]
        output_s = [run_home.fn("g_r.pkl")]
        job.add_task(METHOD, input_s, output_s, n_proc=n_proc, *arg, **kwarg)
    #
    job.to_json()

def read_pkl(pkl_fn_s, group_s):
    dist_bin = []
    _random_distr = []
    _distr_s = {}
    #
    for pkl_fn in pkl_fn_s:
        with pkl_fn.open("rb") as fp:
            X = pickle.load(fp)
        #
        dist_bin.extend(X['dist_bin'])
        _random_distr.extend(X['random'])
        #
        for group_i, group_name_i in group_s.items():
            for group_j, group_name_j in group_s.items():
                distr = []
                for mol_i in group_name_i:
                    for mol_j in group_name_j:
                        mol_pair = (mol_i, mol_j)
                        if mol_pair not in X:
                            continue
                        distr.extend(X[mol_pair])
                if len(distr) == 0:
                    continue
                group_pair = (group_i, group_j)
                if group_pair not in _distr_s:
                    _distr_s[group_pair] = []
                #
                _distr_s[group_pair].append(distr)
    #
    n_bin = [b.shape[0] for b in dist_bin]
    max_bin = max(n_bin)-1
    dist_bin = dist_bin[n_bin.index(max_bin+1)]
    #
    random_distr = np.zeros((len(_random_distr), max_bin), dtype=float)
    for i,x in enumerate(_random_distr):
        random_distr[i,:len(x)] = x
    random_distr = np.mean(random_distr, axis=0)
    nz = (random_distr > 0.)
    #
    distr_s = {}
    for group_pair, Xs in _distr_s.items():
        if len(Xs) == 0:
            continue
        #
        distr_s[group_pair] = []
        for X in Xs:    # iterate traj
            distr = np.zeros_like(random_distr)
            for x in X:
                distr[:len(x)] += x
            distr /= distr.sum()
            distr_s[group_pair].append(distr)
    #
    g_r = {}
    for group_pair, distr in distr_s.items():
        n = len(distr)
        X = np.array(distr)[:,nz] / random_distr[nz][None,:]
        g_r[group_pair] = (np.mean(X, 0), np.std(X, 0)/np.sqrt(n))
    return dist_bin, g_r

def plot(out_home, prefix, group_s, dist_bin, g_r, **kwarg):
    png_fn = out_home.fn(f"g_r.{prefix}.summary.png")
    #
    group_name_s = list(group_s)
    n_group = len(group_name_s)
    fig, axes = plt.subplots(n_group, n_group, figsize=(n_group*3.2, n_group*2.4), sharex=True, sharey=True)
    #
    dist_cntr = 0.5 * (dist_bin[1:] + dist_bin[:-1])
    #
    for group_pair in g_r:
        i = group_name_s.index(group_pair[0])
        j = group_name_s.index(group_pair[1])
        #
        m = g_r[group_pair][0]
        s = g_r[group_pair][1]
        #
        axes[i,j].plot(dist_cntr, m, 'r-', linewidth=1.5)
        axes[i,j].fill_between(dist_cntr, m-s, m+s, color='red', alpha=0.2)
        axes[i,j].plot((0, dist_bin[-1]), (1., 1.), 'k--')
    #
    xlim = kwarg.get("xlim")
    ylim = kwarg.get("ylim")
    for ax in axes.flatten():
        if xlim is not None:
            ax.set_xlim((0, xlim))
        else:
            ax.set_xlim((0, dist_bin[-1]))
        if ylim is not None:
            ax.set_ylim((0, ylim))
        else:
            ax.set_ylim(bottom=0.)
    #
    for i,name in enumerate(group_name_s):
        axes[i,0].set_ylabel(name)
        axes[-1,i].set_xlabel(name)
    #
    fig.tight_layout(h_pad=0.5, w_pad=0.5)
    #
    plt.savefig(png_fn.short())
    plt.close("all")

def summarize(job, prefix, pkl_fn_s, group_s, **kwarg):
    out_fn = job.analysis_home.fn(f"g_r.{prefix}.summary.pkl")
    png_fn = job.analysis_home.fn(f"g_r.{prefix}.summary.png")
    status = (out_fn.status() and png_fn.status())
    if status: return
    #
    dist_bin, g_r = read_pkl(pkl_fn_s, group_s)
    plot(job.analysis_home, prefix, group_s, dist_bin, g_r, **kwarg)
    #
    with out_fn.open("wb") as fout:
        pickle.dump({"dist_bin": dist_bin, "g_r": g_r}, fout)

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
