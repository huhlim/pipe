#!/usr/bin/env python

import os
import sys
import glob
import subprocess as sp

from libcasp14 import WORK_HOME, TARBALL_HOME, SERVER_NAME, CODE_s


def evaluate(ref_fn, model_s):
    cmd = []
    cmd.append("casp_eval.py")
    cmd.extend(["-r", ref_fn])
    cmd.append("-m")
    cmd.extend(model_s)
    cmd.append("-lfast")
    cmd.extend(["-j", "32"])
    #
    output = sp.check_output(cmd).decode("utf8").split("\n")[:-2]
    #
    wrt = []
    dat = []
    for line in output:
        if line.startswith("#"):
            wrt.append("%s\n" % line)
        else:
            value = float(line.strip().split()[1])
            dat.append((value, "%s\n" % line))
    dat.sort(key=lambda x: x[0], reverse=True)
    for x in dat:
        wrt.append(x[1])
    wrt.append("#\n")
    return wrt


def main():
    if len(sys.argv) == 1:
        sys.exit("USAGE: %s [ID]\n" % __file__)
    #
    id = sys.argv[1]
    ref_fn_s = glob.glob("%s/native/%s*.pdb" % (WORK_HOME, id))
    if len(ref_fn_s) == 0:
        sys.exit("No native structure for %s\n" % id)
    ref_fn_s.sort()
    #
    if id.startswith("T"):
        model_s = glob.glob("%s/%s/*_TS?" % (TARBALL_HOME, id))
        model_s = [fn for fn in model_s if not fn.split("/")[-1].startswith("server")]
        model_s.extend(glob.glob("%s/%s/%s/trRosetta/build/min.pdb" % (WORK_HOME, SERVER_NAME, id)))
    elif id.startswith("R"):
        model_s = ["%s/refine.init/%s.pdb" % (WORK_HOME, id)]
    #
    for pred_name in CODE_s.keys():
        if pred_name == SERVER_NAME and id.startswith("T"):
            continue
        model_s.extend(glob.glob("%s/%s/%s/final/model_?.pdb" % (WORK_HOME, pred_name, id)))
    #
    for ref_fn in ref_fn_s:
        if id.startswith("T"):
            out_fn = "%s/%s/%s.dat" % (TARBALL_HOME, id, ref_fn.split("/")[-1][:-4])
        else:
            out_fn = "%s/analysis/TR.raw/%s.dat" % (WORK_HOME, ref_fn.split("/")[-1][:-4])
        if os.path.exists(out_fn):
            continue
        wrt = evaluate(ref_fn, model_s)
        with open(out_fn, "wt") as fout:
            fout.writelines(wrt)


if __name__ == "__main__":
    main()
