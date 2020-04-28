#!/usr/bin/env python

import os
import sys
import copy
import path
import numpy as np
from multiprocessing import Pool

from libtrRosetta import *
from run_trRosetta import read_trRosetta, FEATUREs

PARAM_N_MODEL = 16
EXEC = '%s/apps/trRosetta/scripts/trRosetta.py'%(os.getenv("HOME"))

MAX_MEMORY = 32. # GB

def runner(*args):
    cmd = args[0]
    log = system(cmd, redirect=True, stdout=True)
    return log

def get_energy_min(out_fn_s):
    def read_energy(fn):
        with fn.open() as fp:
            for line in fp:
                if not line.startswith("pose"):
                    continue
                return float(line.strip().split()[-1])
    energy_s = []
    for out_fn in out_fn_s:
        energy_s.append(read_energy(out_fn))
    emin_fn = out_fn_s[np.argmin(energy_s)]
    #
    min_fn = path.Path("min.pdb")
    if min_fn.exists():
        min_fn.remove()
    os.symlink(emin_fn.short(), min_fn.short())
    #
    return min_fn

def run(job, n_model=PARAM_N_MODEL):
    build_home = job.run_home.subdir("build", build=True)
    build_home.chdir()
    #
    npz_fn = job.trRosetta_fn
    fa_fn = job.fa_fn0
    #
    cmd_s = [] ; out_fn_s = []
    for i in range(n_model):
        out_fn = build_home.fn("model.%d.pdb"%i)
        out_fn_s.append(out_fn)
        if out_fn.status():
            continue
        cmd = ['python', EXEC, npz_fn.short(), fa_fn.short(), out_fn.short()]
        cmd_s.append(cmd)
    #
    n_run = len(cmd_s)
    if n_run == 0:
        return get_energy_min(out_fn_s)
    #
    npz = read_trRosetta(npz_fn)
    assert npz['dist'].shape[0] == len(job.seq0), \
            "Sequence and ContactMap do NOT match, %d %d"%(len(job.seq0), npz['dist'].shape[0])
    #
    est_memory = (len(job.seq0)*3000 + 80000)/1024./1024.
    n_iter = int(est_memory * n_run / MAX_MEMORY) + 1
    n_run = int(np.ceil(n_run/n_iter))
    #
    n_proc = min(PARAM_N_PROC, n_run)
    proc = Pool(n_proc)
    log_s = proc.map(runner, cmd_s)
    proc.close()
    #
    with build_home.fn("build.log").open('wt') as fout:
        for log in log_s:
            fout.write(log)
    #
    min_fn = get_energy_min(out_fn_s)

    job.run_home.chdir()
    return min_fn
    
