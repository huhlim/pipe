#!/usr/bin/env python

import os
import sys
import path
import numpy as np

from libtrRosetta import *

METHODs = ["trRosetta"]

OPTIONs = {}
OPTIONs["trRosetta"] = {}
OPTIONs["trRosetta"]["MAX_SEQUENCE"] = 100000
OPTIONs["trRosetta"][
    "HHblits"
] = "-mact 0.35 -maxfilt 1000000 -neffmax 20 -nodiff -realign_max 1000000 -maxseq 1000000 -maxmem 16 -n 4 -v 0".split()
OPTIONs["trRosetta"]["HHblits_database"] = "/feig/s1/huhlim/db/hhsuite/uc30/current/uc30"
OPTIONs["trRosetta"]["HMMbuild"] = "--hand --amino --informat=a2m".split()
OPTIONs["trRosetta"]["HMMsearch"] = "-T 27".split()
OPTIONs["trRosetta"]["HMMsearch_database"] = "/feig/s1/huhlim/db/uniref/uniref100.fasta"


def filter_msa(in_msa, cov, id_cut=90, out_msa=None, max_seq=None):
    if out_msa is None:
        out_msa = path.Path("%s.cov%d" % (in_msa, cov))
        #
    cmd = []
    cmd.append("hhfilter")
    cmd.extend(["-v", "0"])
    cmd.extend(["-i", in_msa.short()])
    cmd.extend(["-o", out_msa.short()])
    cmd.extend(["-id", "%d" % id_cut])
    if cov > 0:
        cmd.extend(["-cov", "%d" % cov])
    cmd.extend(["-maxseq", "1000000"])
    if not out_msa.status():
        system(cmd)
    #
    n_seq = 0
    with out_msa.open() as fp:
        for line in fp:
            if line.startswith(">"):
                n_seq += 1
    if max_seq is None:
        return n_seq
    else:
        if n_seq < max_seq:
            return n_seq
        else:
            id_cut -= 10
            out_msa.remove()
            return filter_msa(in_msa, cov, id_cut=id_cut, out_msa=out_msa, max_seq=max_seq)


def sto2a3m(sto):
    name_s = []
    seq_s = []
    with sto.open() as fp:
        n_chunk = -1
        for line in fp:
            if line.startswith("#"):
                continue
            if line.startswith("//"):
                break
            if line.strip() == "":
                i = -1
                continue
            i += 1
            if i == 0:
                n_chunk += 1
            #
            x = line.strip().split()
            name = x[0]
            seq = x[1]
            if n_chunk == 0:
                name_s.append(name)
                seq_s.append("")
            seq_s[i] += seq

    wrt = []
    for i, name in enumerate(name_s):
        wrt.append(">%s\n" % name)
        wrt.append(seq_s[i].replace(".", "") + "\n")
    return wrt


def run_trRosetta_iterative_HHblits_search(title, fa_fn, n_proc=PARAM_N_PROC):
    has_enough_sequence = False
    #
    eCutoff_s = ["1e-%02d" % e for e in [80, 70, 60, 50, 40, 30, 20, 10, 8, 6, 4]]
    in_msa = fa_fn
    for eCutoff in eCutoff_s:
        out_msa = path.Path("%s.hhblits.%s.a3m" % (title, eCutoff))
        #
        cmd = []
        cmd.append("hhblits")
        cmd.extend(["-i", in_msa.short()])
        cmd.extend(["-oa3m", out_msa.short()])
        cmd.extend(["-o", "/dev/null"])
        cmd.extend(["-cpu", "%d" % n_proc])
        cmd.extend(["-e", eCutoff])
        cmd.extend(["-d", OPTIONs["trRosetta"]["HHblits_database"]])
        cmd.extend(OPTIONs["trRosetta"]["HHblits"])
        if not out_msa.status():
            system(cmd)
        if not out_msa.status():
            sys.exit("Error: failed to generate %s\n" % out_msa)
        #
        n_cov75 = filter_msa(out_msa, 75)
        n_cov50 = filter_msa(out_msa, 50)
        if n_cov75 >= 2000 or n_cov50 >= 5000:
            has_enough_sequence = True
            break
        #
        in_msa = out_msa

    return has_enough_sequence, out_msa


def run_trRosetta_HMMER_search(title, in_msa, n_proc=PARAM_N_PROC):
    hmm_fn = path.Path("%s.hmmer.hmm" % title)
    cmd = []
    cmd.append("hmmbuild")
    cmd.extend(OPTIONs["trRosetta"]["HMMbuild"])
    cmd.append(hmm_fn.short())
    cmd.append(in_msa.short())
    if not hmm_fn.status():
        system(cmd)
    #
    out_sto = path.Path("%s.hmmer.sto" % title)
    cmd = []
    cmd.append("hmmsearch")
    cmd.extend(["-o", "%s.hmmer.log" % title])
    cmd.extend(["-A", out_sto.short()])
    cmd.extend(OPTIONs["trRosetta"]["HMMsearch"])
    cmd.extend(["--cpu", "%d" % n_proc])
    cmd.append(hmm_fn.short())
    cmd.append(OPTIONs["trRosetta"]["HMMsearch_database"])
    if not out_sto.status():
        system(cmd)
    if not out_sto.status():
        sys.exit("Error: failed to generate %s\n" % out_sto)
    #
    out_a3m = path.Path("%s.hmmer.a3m" % title)
    with out_a3m.open("wt") as fout:
        fout.writelines(sto2a3m(out_sto))
    return out_a3m


def run_trRosetta_search(title, fa_fn, n_proc=PARAM_N_PROC):
    final_msa = path.Path("%s.a3m" % title)
    if final_msa.status():
        return final_msa
    #
    has_enough_sequence, out_hhblits_msa = run_trRosetta_iterative_HHblits_search(
        title, fa_fn, n_proc=n_proc
    )
    if has_enough_sequence:
        out_hmmer_msa = None
    else:
        out_hmmer_msa = run_trRosetta_HMMER_search(title, out_hhblits_msa, n_proc=n_proc)
    #
    merged_msa = path.Path("%s.merged.a3m" % title)
    if not merged_msa.status():
        with merged_msa.open("wt") as fout:
            with out_hhblits_msa.open() as fp:
                fout.write(fp.read())
            if out_hmmer_msa is not None:
                with out_hmmer_msa.open() as fp:
                    fout.write(fp.read())
    #
    if not final_msa.status():
        filter_msa(merged_msa, -1, out_msa=final_msa, max_seq=OPTIONs["trRosetta"]["MAX_SEQUENCE"])
    return final_msa


def run(fa_fn, n_proc=PARAM_N_PROC):
    EXECs = {}
    EXECs["trRosetta"] = run_trRosetta_search
    #
    work_home = fa_fn.dirname()
    work_home.chdir()
    #
    title = fa_fn.name()
    #
    msa_s = {}
    for method in METHODs:
        run_home = work_home.subdir(method, build=True)
        run_home.chdir()
        #
        msa = EXECs[method](title, fa_fn)
        msa_s[method] = msa
    return msa_s


def main():
    in_fa = path.Path(sys.argv[1])
    run(in_fa)


if __name__ == "__main__":
    main()
