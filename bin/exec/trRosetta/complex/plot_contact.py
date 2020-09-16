#!/usr/bin/env python

import os
import sys
import path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from importlib import import_module

font={"size": 14}
matplotlib.rc("font", **font)
COLORs = ['red', 'blue', 'green']

D_MAX = 20.0
D_CNT =  8.0
d2c = lambda d: 1. - 1./(1.+np.exp(-2*(d-D_CNT-2)))

def read_pdb(pdb_fn):
    pdb = [] ; resNo_s = [] ; chain_s = []
    with pdb_fn.open() as fp:
        for line in fp:
            if not line.startswith("ATOM"): continue
            #
            resName = line[17:20]
            atmName = line[12:16].strip()
            #
            if resName == 'GLY' and atmName != 'CA':
                continue
            elif resName != 'GLY' and atmName != 'CB':
                continue
            #
            chain = line[21]
            resNo = line[22:26]
            r = np.array([line[30:38], line[38:46], line[46:54]], dtype=float)
            pdb.append(r)
            resNo_s.append(resNo)
            chain_s.append(chain)
    return np.array(pdb), np.array(resNo_s, dtype=int), chain_s

def get_contact_from_PDB(l_seq, pdb_fn, l_fa):
    pdb, resNo_s0, chain_s = read_pdb(pdb_fn)
    #
    k = -1 ; chain_prev = None ; resNo0 = 0
    resNo_s = []
    for resNo, chain in zip(resNo_s0, chain_s):
        if chain != chain_prev:
            k += 1
            chain_prev = chain
            if k > 0:
                resNo0 += l_fa[k-1]
        resNo_s.append(resNo + resNo0)
    #
    status = np.ones(l_seq, dtype=bool)
    R = np.zeros((l_seq, 3), dtype=float)
    k = -1
    for resNo in range(l_seq):
        if np.isin(resNo+1, resNo_s):
            k += 1 ; R[resNo] = pdb[k]
        else:
            status[resNo] = False
    #
    dR = R[:,None] - R[None,:]
    dij = np.sqrt(np.sum(dR**2, axis=-1))
    #
    diso = ~status
    dij[diso] = D_MAX
    dij[:,diso] = D_MAX
    cij = d2c(dij)
    return diso, cij

def run(npz_fn, pdb_fn, png_fn, l_fa):
    npz = np.load(npz_fn.short())
    contact = np.sum(npz['dist'][:,:,1:13], axis=-1)
    l_seq = contact.shape[0]
    #
    diso, pdb = get_contact_from_PDB(l_seq, pdb_fn, l_fa)
    #
    xylim = [0.5, l_seq+0.5]
    args = {"vmin": 0.1, "vmax": 0.25, "cmap": plt.get_cmap("hot_r"), \
            'origin': "lower", "extent": xylim+xylim}
    #
    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    #
    dij = np.zeros((l_seq, l_seq), dtype=float)
    pred = np.triu_indices(l_seq)
    dij[pred] = contact[pred]
    model = np.tril_indices(l_seq)
    dij[model] = pdb[model]
    #
    #cmap0 = ax.imshow(contact, **args)
    cmap0 = ax.imshow(dij, **args)
    ax.plot(xylim, xylim, 'k-')
    ax.grid(True, linestyle='--', color='grey')
    #
    xy = 0.5
    for k,l in enumerate(l_fa):
        color = COLORs[k]
        prot = patches.Rectangle((xy, xy), l, l, linewidth=2, edgecolor=color, facecolor='none')
        ax.add_patch(prot)
        xy += l
    #
    fig.tight_layout()
    plt.savefig(png_fn)
    plt.close("all")

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
        sys.exit("USAGE: %s [NPZ] [PDB] [FAs]\n"%__file__)
    npz_fn = path.Path(sys.argv[1])
    pdb_fn = path.Path(sys.argv[2])
    fa_s = [path.Path(fn) for fn in sys.argv[3:]]
    l_fa = [get_l_seq(fa) for fa in fa_s]
    #
    png_fn = '%s.with_pdb.png'%npz_fn.name()
    #
    run(npz_fn, pdb_fn, png_fn, l_fa)

if __name__ == '__main__':
    main()
