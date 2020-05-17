#!/usr/bin/env python

import os
import sys
import json
import argparse
from importlib import import_module

import warnings
warnings.filterwarnings("ignore")

WORK_HOME = os.getenv("PREFMD_HOME")
assert WORK_HOME is not None
sys.path.insert(0, '%s/bin'%WORK_HOME)

import path
from libcommon import *
from libpdb import Sequence, PDB

def run(arg):
    model_summary = path.Path("model_s")
    if model_summary.status(): return
    #
    model_NA = path.Path("hybrid/DONE")
    model_list = path.Path("hybrid/init_s")
    if not model_list.status() and not model_NA.status():
        cmd = ['%s/prep_hybrid.py'%EXEC_HOME, arg.output_prefix, arg.init_pdb.short()]
        system(cmd)
    #
    if model_NA.status():   # Not enough template
        with model_summary.open("wt") as fout:
            fout.write("#")
        return
    #
    sel_s = path.Path.glob("hybrid/iter_*/sel.out")
    if len(sel_s) < 10:
        cmd = ['%s/run_hybrid.py'%EXEC_HOME, model_list.short(), '%d'%arg.n_proc]
        system(cmd)
    #
    sel_s = path.Path.glob("hybrid/iter_*/sel.out")
    if len(sel_s) == 0:
        with model_summary.open("wt") as fout:
            fout.write("#")
        return
    #
    final = sorted(sel_s, key=lambda x: int(x.split("/")[-2].split("_")[-1]))[-1]
    model_home = path.Dir("model", build=True)
    model_home.chdir()
    #
    cmd = ['%s/main/bin/extract_pdbs.linuxgccrelease'%os.environ['ROSETTA_HOME']]
    cmd.extend(['-in:file:silent', final.short()])
    system(cmd, stdout=True)
    #
    fn_s = sorted(path.Path.glob("iter*.*.pdb"), key=lambda x: int(x.split(".")[-2]))
    model_s = []
    for i,fn in enumerate(fn_s):
        out_fn = path.Path("model.%d.pdb"%i)
        with out_fn.open("wt") as fout:
            system(['%s/match_resNo.py'%EXEC_HOME, arg.init_pdb.short(), fn.short()], outfile=fout)
        model_s.append(out_fn)
    #
    with model_summary.open("wt") as fout:
        for model in model_s:
            fout.write("%s\n"%model)

def main():
    arg = argparse.ArgumentParser(prog='hybrid')
    arg.add_argument(dest='output_prefix')
    arg.add_argument(dest='init_pdb')
    arg.add_argument('-j', '--cpu', dest='n_proc', type=int)
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    arg.init_pdb = path.Path(arg.init_pdb)
    #
    run(arg)

if __name__=='__main__':
    main()
