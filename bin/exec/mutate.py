#!/usr/bin/env python

import os
import sys
import json
import argparse

WORK_HOME = os.getenv("PREFMD_HOME")
assert WORK_HOME is not None
sys.path.insert(0, '%s/bin'%WORK_HOME)

from libcommon import *

import seqName

def write_fa(pdb_fn, mutation_s):
    residue_s = []
    seq = []
    with open(pdb_fn) as fp:
        for line in fp:
            if not (line.startswith("ATOM") or line.startswith("HETA")):
                if line.startswith("END"):
                    break
                continue
            if line[12:16].strip() != 'CA':
                continue
            if line[16] not in [' ', 'A']:
                continue
            resNo = line[22:26].strip()
            chain_id = line[21].strip()
            residue_s.append((resNo, chain_id))
            seq.append(seqName.to_one_letter(line[17:20].strip()))
    #
    for mut in mutation_s:
        try:
            key,aa = mut.split(":")
            key = key.split(".")
            resNo = key[0]
            if len(key) == 1:
                chain = None
            elif len(key) == 2:
                chain = key[1]
            else:
                raise
        except:
            sys.exit("Incorrect mutation argument %s\n"%mut)
        #
        if (resNo, chain) in residue_s:
            index = residue_s.index((resNo, chain))
            seq[index] = aa
        elif (chain is None) and (resNo, "") in residue_s:
            index = residue_s.index((resNo, ""))
            seq[index] = aa
        else:
            sys.exit("Failed to find %s\n"%mut)

    return ''.join(seq)

def main():
    arg = argparse.ArgumentParser(prog='mutate')
    arg.add_argument('-i', '--input', dest='input_pdb', required=True)
    arg.add_argument('-o', '--output', dest='output_pdb', required=True)
    arg.add_argument('-m', '--mutate', dest='mutation_s', nargs='+', default=[])
    # --mutate 12:F or --mutate 12.A:F
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    #
    cmd = []
    cmd.append("scwrl4")
    cmd.extend(['-i', arg.input_pdb])
    cmd.extend(['-o', arg.output_pdb])
    #
    if len(arg.mutation_s) > 0:
        seq = write_fa(arg.input_pdb, arg.mutation_s)
        fa_fn = path.Path(path.prefix(arg.input_pdb) + '.fa')
        with fa_fn.open("wt") as fout:
            fout.write(seq)
        cmd.extend(['-s', fa_fn.short()])
    else:
        fa_fn = None
    #
    system(cmd, verbose=False, stdout=True, errfile='/dev/null')
    #
    if fa_fn is not None:
        fa_fn.remove()

if __name__ == '__main__':
    main()
