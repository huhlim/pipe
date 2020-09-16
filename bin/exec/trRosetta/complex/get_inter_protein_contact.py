#!/usr/bin/env python

import os
import sys
import path
import numpy as np
from itertools import combinations

PROB_CUTOFF = 0.10
DISTs = np.arange(36)*0.5+2.25

def run(npz_fn, txt_fn, l_fa):
    npz = np.load(npz_fn.short())
    contact0 = npz['dist'][:,:,1:]
    contact = np.sum(npz['dist'][:,:,1:13], axis=-1)
    #
    i_res = np.cumsum(l_fa)
    n_prot = len(l_fa)-1
    #
    wrt = []
    for i,j in combinations(range(n_prot), 2):
        cnt_ij0 = contact0[i_res[i]:i_res[i+1], i_res[j]:i_res[j+1]]
        cnt_ij = contact[i_res[i]:i_res[i+1], i_res[j]:i_res[j+1]]
        #
        pair_s = np.where(cnt_ij > PROB_CUTOFF)
        dij_s = DISTs[np.argmax(cnt_ij0[pair_s], axis=-1)]
        prob_s = cnt_ij[pair_s]
        index = np.argsort(prob_s)[::-1]
        for k in index:
            res_i = pair_s[0][k]+1
            res_j = pair_s[1][k]+1
            dij = dij_s[k]
            prob = prob_s[k]
            wrt.append("%d-%d %4d %4d %6.4f %5.2f\n"%(i+1, j+1, res_i, res_j, prob, dij))
    with open(txt_fn, 'wt') as fout:
        fout.writelines(wrt)

def get_l_seq(fa_fn):
    seq = []
    with fa_fn.open() as fp:
        for line in fp:
            if not line.startswith(">"):
                seq.append(line.strip())
    seq = ''.join(seq)
    l_seq = len(seq)
    return l_seq

def main():
    if len(sys.argv) < 4:
        sys.exit("USAGE: %s [NPZ] [FAs]\n"%__file__)
    npz_fn = path.Path(sys.argv[1])
    fa_s = [path.Path(fn) for fn in sys.argv[2:]]
    l_fa = [0] + [get_l_seq(fa) for fa in fa_s]
    #
    txt_fn = '%s.txt'%npz_fn.name()
    #
    run(npz_fn, txt_fn, l_fa)

if __name__ == '__main__':
    main()
