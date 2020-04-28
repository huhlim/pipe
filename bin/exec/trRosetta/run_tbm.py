#!/usr/bin/env python

import os
import sys
import numpy as np

import path
from libcommon import *

METHODs = ['blast', 'hhpred', 'hmmer']

EXECs = {}
EXECs['blast']  = 'run_blast.py'
EXECs['hhpred'] = 'run_hhpred.py'
EXECs['hmmer']  = 'run_hmmer.py'

OPTIONs = {}
OPTIONs['blast']  = ['--allow_fix_build']
OPTIONs['hhpred'] = ['--allow_fix_build']
OPTIONs['hmmer']  = ['--allow_fix_build']

def run(fa_fn, n_build=10, n_proc=8, exclude=TBM_EXCLUDE):
    work_home = fa_fn.dirname()
    work_home.chdir()
    #
    title = fa_fn.name()
    #
    summary_s = {}
    for method in METHODs:
        run_home = work_home.subdir(method, build=True)
        run_home.chdir()
        #
        out_fn = run_home.fn("%s.templ_s.summary"%title)
        #
        cmd = []
        cmd.append(EXECs[method])
        cmd.append(fa_fn.short())
        cmd.extend(['--protocol', 'build'])
        cmd.extend(['--n_templ', '%d'%n_build])
        cmd.extend(['--cpu', '%d'%n_proc])
        if exclude is not None:
            cmd.extend(['--exclude', exclude])
        cmd.extend(OPTIONs[method])
        #
        if not out_fn.status():
            system(cmd)
        #
        summary_s[method] = read_summary(out_fn)
    return summary_s

def read_tbm(title, work_home):
    summary_s = {}
    for method in METHODs:
        run_home = work_home.subdir(method)
        out_fn = run_home.fn("%s.templ_s.summary"%title)
        if out_fn.status():
            summary_s[method] = read_summary(out_fn)
        else:
            summary_s[method] = []
    return summary_s

def read_summary(fn):
    summary = []
    with fn.open() as fp:
        for line in fp:
            x = line.strip().split()
            summary.append((float(x[0]), float(x[1]), float(x[2]), x[3], x[4]))
    return summary

def main():
    in_fa = path.Path(sys.argv[1])
    run(in_fa)

if __name__ == '__main__':
    main()
