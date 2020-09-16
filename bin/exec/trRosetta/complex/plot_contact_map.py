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

def run(npz_fn, png_fn, l_fa):
    npz = np.load(npz_fn.short())
    contact = np.sum(npz['dist'][:,:,1:13], axis=-1)
    l_seq = contact.shape[0]
    #
    xylim = [0.5, l_seq+0.5]
    args = {"vmin": 0.1, "vmax": 0.25, "cmap": plt.get_cmap("hot_r"), \
            'origin': "lower", "extent": xylim+xylim}
    #
    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    #
    cmap0 = ax.imshow(contact, **args)
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
        sys.exit("USAGE: %s [NPZ] [FAs]\n"%__file__)
    npz_fn = path.Path(sys.argv[1])
    fa_s = [path.Path(fn) for fn in sys.argv[2:]]
    l_fa = [get_l_seq(fa) for fa in fa_s]
    #
    png_fn = '%s.png'%npz_fn.name()
    #
    run(npz_fn, png_fn, l_fa)

if __name__ == '__main__':
    main()
