#!/usr/bin/env python

import os
import sys
import numpy as np

L_SPACER = 20
def read_fa(fn):
    seq = []
    with open(fn) as fp:
        for line in fp:
            if not line.startswith(">"):
                seq.append(line.strip())
    seq = ''.join(seq)
    l_seq = len(seq)
    return l_seq

def main():
    if len(sys.argv) < 5:
        sys.exit("USAGE: %s [in_NPZ] [out_NPZ] [FAs]\n"%__file__)
    #
    in_npz_fn = sys.argv[1]
    out_npz_fn = sys.argv[2]
    fa_s = sys.argv[3:]
    #
    in_npz = np.load(in_npz_fn)
    mask = np.ones(in_npz['dist'].shape[0], dtype=bool)
    #
    res0 = 0
    for fa in fa_s[:-1]:
        l_seq = read_fa(fa)
        mask[l_seq+res0:l_seq+res0+L_SPACER] = False
        res0 += l_seq + L_SPACER
    MASK = np.ix_(mask, mask)
    #
    out = {}
    for feature in in_npz:
        out[feature] = in_npz[feature][MASK]
    np.savez(out_npz_fn, **out)

if __name__ == '__main__':
    main()
