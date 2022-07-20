#!/usr/bin/env python

import os
import sys
import path
import argparse
import numpy as np
import mdtraj
from multiprocessing import Pool

from libseq import Sequence
from libtm import TM_score

from run_trRosetta import FEATUREs, read_trRosetta
from tbm_to_contact import get_features, get_distr
from build_model import get_energy_min, runner

PARAM_N_PROC = 16
PARAM_N_MODEL = 16
MAX_MEMORY = 32.0  # GB
EXEC = "%s/apps/trRosetta/scripts/trRosetta.py" % (os.getenv("HOME"))


def extract_contacts(weight_s, distr_s):
    out = {}
    for feature in FEATUREs:
        out[feature] = np.zeros_like(distr_s[0][feature])
        l_seq = out[feature][0].shape[0]
        weight_sum = np.zeros((l_seq, l_seq), dtype=np.float32)
        #
        for weight, distr in zip(weight_s, distr_s):
            w = np.ix_(weight[1], weight[1])
            weight_sum[w] += weight[0]
            out[feature] += distr[feature] * weight[0]
        out[feature] /= weight_sum[:, :, None]
        zero = weight_sum == 0.0
        out[feature][zero] = 0.0
    return out, weight_sum


def hybrid_contacts(npz_s):
    hybrid = {}
    for feature in FEATUREs:
        hybrid[feature] = np.zeros_like(npz_s[0][1][feature])
        l_seq = hybrid[feature][0].shape[0]
        weight_sum = np.zeros((l_seq, l_seq), dtype=np.float32)
        #
        for w, npz in npz_s:
            hybrid[feature] += w[:, :, None] * npz[feature]
            weight_sum += w
        hybrid[feature] /= weight_sum[:, :, None]
    return hybrid


def model_to_contact(pdb_fn_s, weight_s):
    distr_s = []
    for weight, pdb_fn in zip(weight_s, pdb_fn_s):
        pdb = mdtraj.load(pdb_fn.short())
        pdb = pdb.atom_slice(pdb.top.select("name CA or name CB or name N or name O"))
        #
        n_residue = pdb.top.n_residues
        is_gly, feature_s = get_features(pdb)
        distr_s.append(get_distr(n_residue, feature_s, weight[1]))
    #
    out, weight_sum = extract_contacts(weight_s, distr_s)
    return out, weight_sum


def prep(arg):
    fa_fn = path.Path(arg.fa_fn)
    pdb_fn_s = [path.Path(fn) for fn in arg.pdb_fn_s]
    #
    seq = Sequence.read(fa_fn.short())
    #
    wrt = []
    wrt.append("# SEQUENCE %s\n" % fa_fn)
    wrt.append("%s\n" % seq.fasta)
    #
    tm_s = []
    for pdb_fn in pdb_fn_s:
        tm = TM_score(pdb_fn, arg.init_fn).tm
        tm_s.append((tm**2, pdb_fn))
    tm_s.sort(key=lambda x: x[0], reverse=True)
    tm_sum = sum([tm[0] for tm in tm_s])
    #
    for tm, pdb_fn in tm_s:
        wrt.append("# STRUCTURE %6.4f %s\n" % (tm / tm_sum, pdb_fn))
        wrt.append("%s\n" % seq.fasta)
    #
    with open("hybrid.in", "wt") as fout:
        fout.writelines(wrt)


def read_input(fn):
    sequence = None
    weight_s = []
    pdb_fn_s = []
    with fn.open() as fp:
        for line in fp:
            if line.startswith("#"):
                x = line.strip().split()[1:]
                if x[0] == "SEQUENCE":
                    fa_fn = path.Path(x[1])
                else:
                    weight = [float(x[1])]
                    weight_s.append(weight)
                    pdb_fn_s.append(path.Path(x[2]))
            else:
                if sequence is None:
                    sequence = np.array(list(line.strip()), dtype="<U1") != "-"
                else:
                    mask = np.array(list(line.strip()), dtype="<U1") != "-"
                    weight.append(mask)
    #
    weight_sum = sum([tm[0] for tm in weight_s])
    for i, w in enumerate(weight_s):
        weight_s[i][0] = w[0] / weight_sum
    #
    return fa_fn, sequence, weight_s, pdb_fn_s


def build_model(npz_fn, fa_fn, l_seq, n_model=PARAM_N_MODEL):
    cmd_s = []
    out_fn_s = []
    for i in range(n_model):
        out_fn = path.Path("model.%d.pdb" % i)
        out_fn_s.append(out_fn)
        if out_fn.status():
            continue
        cmd = ["python", EXEC, npz_fn.short(), fa_fn.short(), out_fn.short()]
        cmd_s.append(cmd)
    #
    n_run = len(cmd_s)
    if n_run == 0:
        return get_energy_min(out_fn_s)
    #
    npz = read_trRosetta(npz_fn)
    assert npz["dist"].shape[0] == l_seq, "Sequence and ContactMap do NOT match, %d %d" % (
        l_seq,
        npz["dist"].shape[0],
    )
    #
    est_memory = (l_seq * 3000 + 80000) / 1024.0 / 1024.0
    n_iter = int(est_memory * n_run / MAX_MEMORY) + 1
    n_run = int(np.ceil(n_run / n_iter))
    #
    n_proc = min(PARAM_N_PROC, n_run)
    proc = Pool(n_proc)
    log_s = proc.map(runner, cmd_s)
    proc.close()
    #
    with path.Path("build.log").open("wt") as fout:
        for log in log_s:
            fout.write(log)
    #
    min_fn = get_energy_min(out_fn_s)


def run(arg):
    arg.input_fn = path.Path(arg.input_fn)
    arg.npz_fn = path.Path(arg.npz_fn)
    #
    fa_fn, sequence, weight_s, pdb_fn_s = read_input(arg.input_fn)
    l_seq = len(sequence)
    #
    npz0 = read_trRosetta(arg.npz_fn)
    npz0_weight = np.zeros((l_seq, l_seq), dtype=np.float32)
    npz0_weight[np.ix_(sequence, sequence)] += arg.npz_weight
    #
    model_npz, weight_sum = model_to_contact(pdb_fn_s, weight_s)
    weight_sum *= 1.0 - arg.npz_weight
    #
    npz_s = []
    npz_s.append((npz0_weight, npz0))
    npz_s.append((weight_sum, model_npz))
    hybrid = hybrid_contacts(npz_s)
    #
    npz_fn = path.Path("hybrid.npz")
    np.savez(npz_fn.short(), **hybrid)
    #
    build_model(npz_fn, fa_fn, l_seq)


def main():
    arg = argparse.ArgumentParser(prog="hybrid_models")
    sub = arg.add_subparsers(dest="method")
    #
    arg_prep = sub.add_parser("prep")
    arg_prep.add_argument("--init", dest="init_fn", required=True)
    arg_prep.add_argument("--seq", dest="fa_fn", required=True)
    arg_prep.add_argument("--model", dest="pdb_fn_s", nargs="*", default=[])
    #
    arg_run = sub.add_parser("run")
    arg_run.add_argument("--init", dest="init_fn", required=True)
    arg_run.add_argument("--input", dest="input_fn", required=True)
    arg_run.add_argument("--npz", dest="npz_fn", required=True)
    arg_run.add_argument("--weight", dest="npz_weight", type=float, default=0.5)

    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args()
    #
    arg.init_fn = path.Path(arg.init_fn)
    #
    if arg.method == "prep":
        prep(arg)
    elif arg.method == "run":
        run(arg)


if __name__ == "__main__":
    main()
