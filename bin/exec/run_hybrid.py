#!/usr/bin/env python

import os
import sys
import numpy as np
import subprocess as sp

WORK_HOME = os.getenv("PIPE_HOME")
assert WORK_HOME is not None
sys.path.insert(0, "%s/bin" % WORK_HOME)

import path
from libcommon import *

EXEC_HOME = path.Path(__file__).dirname()
EXEC_RESTRAINT = EXEC_HOME.fn("generate_harmonic_cst.py")

ROSETTA_HOME = os.getenv("ROSETTA_HOME")
SCRIPT_HOME = "%s/main/source/scripts/python/public/iterative_hybridize" % ROSETTA_HOME
PSIPRED_EXEC = "%s/apps/psipred/current/run_psipred" % os.getenv("HOME")

VERBOSE = True


def prep_init_models(init_s):
    prep_s = []
    init = path.Path("init.pdb")
    prep_s.append(init)
    if not init.status():
        with init.open("wt") as fout:
            system(
                "convpdb.pl -out generic -renumber 1 %s" % init_s[0].short(),
                outfile=fout,
                verbose=False,
            )
    #
    disu_s = []
    with init.open() as fp:
        for line in fp:
            if line.startswith("SSBOND"):
                disu_s.append("%s %s\n" % (line[17:21], line[31:35]))
    with open("disulf.def", "wt") as fout:
        if len(disu_s) > 0:
            fout.writelines(disu_s)
        else:
            fout.write("1 1\n")
    #
    for i, model in enumerate(init_s[1:]):
        prep = path.Path("model.%d.pdb" % (i + 1))
        if not prep.status():
            with prep.open("wt") as fout:
                system(["%s/match_resNo.py" % EXEC_HOME, init.short(), model.short()], outfile=fout)
        prep_s.append(prep)
    return prep_s


def build_frag(fa_fn, n_proc):
    name = fa_fn.name()
    #
    if not os.path.exists("%s.200.3mers" % name) or not os.path.exists("%s.200.9mers" % name):
        ss2_fn = path.Path("%s.ss2" % name)
        if not ss2_fn.status():
            cmd = [PSIPRED_EXEC, fa_fn.short()]
            system(cmd, errfile="/dev/null")
        #
        # cmd = ['python2', "%s/tools/fragment_tools/make_fragments.py"%ROSETTA_HOME]
        cmd = ["python", "%s/tools/fragment_tools/make_fragments.py" % ROSETTA_HOME]
        cmd.extend(["-cpus", "%d" % n_proc])
        cmd.extend(["-psipredfile", ss2_fn.path()])
        cmd.append(fa_fn.short())
        system(cmd, errfile="/dev/null")
    if not os.path.exists("%s.200.3mers" % name) or not os.path.exists("%s.200.9mers" % name):
        return False
    if not os.path.exists("t000_.3mers"):
        os.symlink("%s.200.3mers" % name, "t000_.3mers")
    if not os.path.exists("t000_.9mers"):
        os.symlink("%s.200.9mers" % name, "t000_.9mers")
    return True


def prep_run(init_s):
    if not os.path.exists("model.out"):
        cmd = []
        cmd.append("%s/main/source/bin/combine_silent.linuxgccrelease" % ROSETTA_HOME)
        cmd.append("-in:file:s")
        cmd.extend([fn.short() for fn in init_s])
        cmd.extend(["-out:file:silent", "model.out"])
        cmd.extend(["-out:file:silent_struct_type", "binary"])
        cmd.extend(["-ignore_zero_occupancy", "false"])
        system(cmd, verbose=True)
    #
    if not os.path.exists("ref.out"):
        cmd = []
        cmd.append("%s/main/source/bin/iterhybrid_selector.linuxgccrelease" % ROSETTA_HOME)
        cmd.extend(["-in:file:silent", "model.out"])
        cmd.extend(["-in:file:template_pdb", "init.pdb"])
        cmd.extend(["-cm:similarity_cut", "0.2"])
        cmd.extend(["-out:file:silent ref.out"])
        cmd.extend(["-out:nstruct", "%d" % len(init_s)])
        cmd.extend(["-out:prefix", "iter0"])
        cmd.extend(["-score:weights", "ref2015_cart"])
        cmd.extend(["-silent_read_through_errors"])
        cmd.extend(["-in:file:silent_struct_type", "binary"])
        cmd.extend(["-out:file:silent_struct_type", "binary"])
        cmd.extend(["-mute", "core", "basic"])
        system(cmd, verbose=True)
    #
    if not os.path.exists("cen.cst"):
        cmd = []
        cmd.append("%s/generate_harmonic_cst.py" % EXEC_HOME)
        cmd.extend([fn.short() for fn in init_s])
        with open("cen.cst", "wt") as fout:
            system(cmd, verbose=True, outfile=fout)
    if not os.path.exists("fa.cst"):
        os.symlink("cen.cst", "fa.cst")


def main():
    init_fn_s = path.Path(sys.argv[1])
    if len(sys.argv) > 2:
        n_proc = int(sys.argv[2])
    else:
        n_proc = 20
    #
    hybrid_home = init_fn_s.dirname()
    hybrid_home.chdir()
    #
    init_s = []
    with init_fn_s.open() as fp:
        for line in fp:
            init_s.append(path.Path(line.strip()))
    #
    init_s = prep_init_models(init_s)
    #
    fa_fn0 = hybrid_home.fn("input0.fa")
    fa_fn = hybrid_home.fn("input.fa")
    if (not fa_fn0.status()) or (not fa_fn.status()):
        out_s = system(
            ["%s/pdb_seq" % EXEC_HOME, init_s[0].short()], stdout=True, verbose=False
        ).split("\n")[:-1]
        with fa_fn0.open("wt") as fout:
            fout.write(">init\n")
            for line in out_s:
                if not line.startswith(">"):
                    fout.write("%s\n" % line)
        with fa_fn.open("wt") as fout:
            fout.write(">init\n")
            for line in out_s[1:]:
                if line.startswith(">"):
                    fout.write("/\n")
                else:
                    fout.write("%s" % line)
            fout.write("\n")
    #
    if not build_frag(fa_fn0, n_proc):
        sys.exit("Failed to build fragments\n")
    #
    prep_run(init_s)
    #
    with open("NODEFILE", "wt") as fout:
        fout.write("%d/localhost" % n_proc)
    #
    cmd = []
    # cmd.append("python2")
    cmd.append("python")
    cmd.append("%s/IterationMaster.py" % SCRIPT_HOME)
    cmd.extend(["-iha", "60"])
    cmd.extend(["-nodefile", "NODEFILE"])
    cmd.append("-simple")
    cmd.append("-niter")
    cmd.append("10")
    with open("iter_hybrid.log", "wt") as fout:
        system(cmd, outfile=fout, errfile=fout, verbose=True)


if __name__ == "__main__":
    main()
