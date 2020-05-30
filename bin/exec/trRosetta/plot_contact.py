#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from importlib import import_module

font={"size": 14}
matplotlib.rc("font", **font)

D_MAX = 20.0
D_CNT =  8.0
d2c = lambda d: 1. - 1./(1.+np.exp(-2*(d-D_CNT-2)))

def read_pdb(pdb_fn):
    pdb = [] ; resNo_s = []
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
            resNo = line[22:26]
            r = np.array([line[30:38], line[38:46], line[46:54]], dtype=float)
            pdb.append(r)
            resNo_s.append(resNo)
    return np.array(pdb), np.array(resNo_s, dtype=int)

def get_contact_from_PDB(l_seq, pdb_fn):
    pdb, resNo_s = read_pdb(pdb_fn)
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

def run(npz_fn, png_fn0, pdb_fn_s):
    npz = np.load(npz_fn.short())
    contact = np.sum(npz['dist'][:,:,1:13], axis=-1)
    l_seq = contact.shape[0]
    #
    xylim = [0.5, l_seq+0.5]
    args = {"vmin": 0.0, "vmax": 0.25, "cmap": plt.get_cmap("hot_r"), \
            'origin': "lower", "extent": xylim+xylim}
    #
    for pdb_fn in pdb_fn_s:
        if png_fn0 is None:
            png_fn = '%s.png'%(pdb_fn[:-4])
        else:
            png_fn = png_fn0.short()
        #
        diso, pdb = get_contact_from_PDB(l_seq, pdb_fn)
        #
        fig, ax = plt.subplots(figsize=(6.4, 5.4))
        #
        dij = np.zeros((l_seq, l_seq), dtype=float)
        pred = np.triu_indices(l_seq)
        dij[pred] = contact[pred]
        model = np.tril_indices(l_seq)
        dij[model] = pdb[model]
        #
        dij[diso] = 0.0
        dij[:,diso] = 0.0
        #
        score_pred = contact[np.triu_indices(l_seq, k=12)]
        score_pred.sort()
        score_pred = score_pred[::-1]
        score_pred = np.mean(score_pred[:l_seq])
        #
        score_model = pdb[np.triu_indices(l_seq, k=12)]
        score_model.sort()
        score_model = score_model[::-1]
        score_model = np.mean(score_model[:int(l_seq*1.5)])
        #
        cmap0 = ax.imshow(dij, **args)
        ax.plot(xylim, xylim, 'k-')
        ax.grid(True, linestyle='--', color='grey')
        #
        ax.set_xlabel("from Contact Prediction (score=%6.4f)"%score_pred)
        ax.set_ylabel("from PDB (score=%6.4f)"%score_model)
        #
        fig.tight_layout()
        plt.savefig(png_fn)
        plt.close("all")

def main():
    path = import_module("path")
    npz_fn = path.Path(sys.argv[1])
    #
    png_fn = None
    pdb_fn_s = []
    for arg in sys.argv[2:]:
        if arg.endswith(".png"):
            png_fn = path.Path(arg)
        else:
            pdb_fn_s.append(path.Path(arg))
    #
    run(npz_fn, png_fn, pdb_fn_s)

if __name__ == '__main__':
    main()
