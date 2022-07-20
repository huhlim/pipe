#!/usr/bin/env python

import os
import sys
import numpy as np
import subprocess as sp

WORK_HOME = os.getenv("PIPE_HOME")
assert WORK_HOME is not None
sys.path.insert(0, "%s/bin" % WORK_HOME)

from libcommon import *

import path
import libhhsuite

N_PROC = 12
N_TEMPL = 100
N_HYBRID = 10
N_HYBRID_MIN = 2
TM_CUTOFF = 0.6
TM_OFFSET = 0.2

VERBOSE = False
EXCLUDE = []
if TBM_EXCLUDE is not None:
    with open(TBM_EXCLUDE) as fp:
        for line in fp:
            EXCLUDE.append(line.strip())


def calc_tmalign(ref_fn, pdb_s):
    cmd = []
    cmd.append("%s/calc_tmalign" % EXEC_HOME)
    cmd.append("--cpu")
    cmd.append("%d" % N_PROC)
    cmd.append("--ref")
    cmd.append(ref_fn.short())
    cmd.append("--list")
    cmd.append(pdb_s)
    lines = system(cmd, stdout=True, verbose=False).split("\n")[:-1]
    score = [float(line.strip().split()[1]) for line in lines]
    return score


def calc_tmscore(ref_fn, pdb_s):
    cmd = []
    cmd.append("%s/calc_tmscore" % EXEC_HOME)
    cmd.append("--cpu")
    cmd.append("%d" % N_PROC)
    cmd.append("--ref")
    cmd.append(ref_fn.short())
    cmd.append("--list")
    cmd.append(pdb_s)
    lines = system(cmd, stdout=True, verbose=False).split("\n")[:-1]
    score = [float(line.strip().split()[1]) for line in lines]
    return score


def run_hhsearch(homolog_home, id, fa_fn, input_pdb):
    hh_fn = homolog_home.fn("%s.vit.local" % id)
    if not hh_fn.status():
        cmd = ["%s/run_hhpred.py" % EXEC_HOME, fa_fn.short()]
        cmd.extend(["-p", "simple"])
        cmd.extend(["-d", HH_sequence_database])
        print(" ".join(cmd))
        system(cmd, verbose=False, stdout=True, errfile="/dev/null")
    #
    hh_s = libhhsuite.parse_hhr(hh_fn.short())
    #
    tm_fn = homolog_home.fn("tm.dat")
    if not tm_fn.status():
        pdb_fn_s = []
        for hh in hh_s[:N_TEMPL]:
            pdb_id = hh.pdb_id
            if pdb_id[:4] in EXCLUDE:
                continue
            #
            pdb_fn = homolog_home.fn("%s.pdb" % pdb_id)
            if not pdb_fn.status():
                system(
                    ["%s/pdb_get" % EXEC_HOME, "-f", pdb_id],
                    verbose=False,
                    stdout=True,
                    errfile="/dev/null",
                )
            if pdb_fn.status():
                pdb_fn_s.append(pdb_fn)
        with homolog_home.fn("pdb_s").open("wt") as fout:
            for pdb_fn in pdb_fn_s:
                fout.write(pdb_fn.short() + "\n")
        #
        tm_s = calc_tmalign(input_pdb, "pdb_s")
        selected = []
        with tm_fn.open("wt") as fout:
            for pdb_fn, tm in zip(pdb_fn_s, tm_s):
                fout.write("%s  " % pdb_fn.name())
                fout.write(" %6.4f" % tm)
                fout.write("\n")
                if tm > TM_CUTOFF:
                    if pdb_fn.name() not in selected:
                        selected.append(pdb_fn.name())
                pdb_fn.remove()
        with homolog_home.fn("selected").open("wt") as fout:
            for pdb_id in selected:
                fout.write("%s\n" % pdb_id)
    else:
        selected = []
        with homolog_home.fn("selected").open() as fp:
            for line in fp:
                selected.append(line.strip())
    return selected


def run_modeller(homolog_home, id, fa_fn, input_pdb, selected):
    model_fn = homolog_home.fn("model_s")
    if not model_fn.status():
        cmd = []
        cmd.append("%s/run_hhpred.py" % EXEC_HOME)
        cmd.append(fa_fn.short())
        cmd.append("-p")
        cmd.append("build")
        cmd.append("--include")
        cmd.append("selected")
        cmd.append("--force")
        cmd.append("--alt")
        cmd.append("3")
        cmd.append("--cpu")
        cmd.append("%d" % N_PROC)
        cmd.append("-d")
        cmd.append(HH_sequence_database)
        cmd.append("--hhdb")
        cmd.append(HH_pdb70_database)
        if len(selected) > 0:
            print(" ".join(cmd))
            system(cmd, verbose=False, stdout=True, errfile="/dev/null")
        #
        model_s = []
        for pdb_id in selected:
            model = path.Path.glob("%s-%s*.model.pdb" % (id, pdb_id))
            for m in model:
                name = ".".join(m.fname().split(".")[:-2])
                out_fn = homolog_home.fn("%s.pdb" % name)
                if not out_fn.status():
                    cmd = []
                    cmd.append("%s/match_resNo.py" % EXEC_HOME)
                    cmd.append(input_pdb.short())
                    cmd.append(m.short())
                    with out_fn.open("wt") as fout:
                        system(cmd, outfile=fout, verbose=False, stdout=True, errfile="/dev/null")
                if out_fn not in model_s:
                    model_s.append(out_fn)

        with model_fn.open("wt") as fout:
            for model in model_s:
                fout.write("%s\n" % model.short())
    else:
        model_s = []
        with model_fn.open() as fp:
            for line in fp:
                model_s.append(path.Path(line.strip()))
    #
    tm_fn = homolog_home.fn("model.dat")
    if not tm_fn.status():
        tm_s = calc_tmscore(input_pdb, "model_s")
        with tm_fn.open("wt") as fout:
            for pdb_fn, tm in zip(model_s, tm_s):
                fout.write("%6.4f %s\n" % (tm, pdb_fn))
    else:
        tm_s = []
        with tm_fn.open() as fp:
            for line in fp:
                x = line.strip().split()
                tm_s.append(float(x[0]))
    return model_s, tm_s


def select_init(hybrid_home, input_pdb, model_s, tm_s):
    if len(model_s) > 0:
        tm_cutoff = max(max(tm_s) - TM_OFFSET, TM_CUTOFF)
        sorted_index = np.argsort(tm_s)[::-1]
        init_s = []
        init_s.append(input_pdb)
        for i in sorted_index:
            tm = tm_s[i]
            model = model_s[i]
            if tm > tm_cutoff:
                init_s.append(model)
        init_s = init_s[:N_HYBRID]
    else:
        init_s = [input_pdb]
    #
    hybrid_home.chdir()
    if len(init_s) < N_HYBRID_MIN:
        with hybrid_home.fn("DONE").open("wt") as fout:
            fout.write("Failed to detect templates for hybridization\n")
    else:
        with hybrid_home.fn("init_s").open("wt") as fout:
            n_init = len(init_s)
            for i in range(N_HYBRID):
                fout.write("%s\n" % init_s[i % n_init].short())


def use_selected(hybrid_home, input_pdb, selected_home):
    init_s = []
    init_s.append(input_pdb)
    init_s.extend(selected_home.glob("*.pdb"))
    #
    hybrid_home.chdir()
    with hybrid_home.fn("init_s").open("wt") as fout:
        n_init = len(init_s)
        for i in range(N_HYBRID):
            fout.write("%s\n" % init_s[i % n_init].short())


def main():
    id = sys.argv[1]
    input_pdb = path.Path(sys.argv[2])
    #
    selected_home = path.Dir("selected")
    homolog_home = path.Dir("homolog", build=True)
    hybrid_home = path.Dir("hybrid", build=True)
    #
    if not selected_home.status():
        #
        homolog_home.chdir()
        #
        fa_fn = homolog_home.fn("%s.fa" % id)
        if not fa_fn.status():
            with fa_fn.open("wt") as fout:
                system(
                    ["%s/pdb_seq" % EXEC_HOME, input_pdb.short()],
                    outfile=fout,
                    verbose=False,
                    stdout=True,
                    errfile="/dev/null",
                )
        #
        selected = run_hhsearch(homolog_home, id, fa_fn, input_pdb)
        #
        model_s, tm_s = run_modeller(homolog_home, id, fa_fn, input_pdb, selected)
        #
        select_init(hybrid_home, input_pdb, model_s, tm_s)
    else:
        use_selected(hybrid_home, input_pdb, selected_home)


if __name__ == "__main__":
    main()
