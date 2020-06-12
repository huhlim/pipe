#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from importlib import import_module

font={"size": 14}
matplotlib.rc("font", **font)

def run(npz_fn, png_fn):
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
    fig.tight_layout()
    plt.savefig(png_fn)
    plt.close("all")

def main():
    path = import_module("path")
    npz_fn = path.Path(sys.argv[1])
    #
    png_fn = '%s.png'%npz_fn.name()
    #
    run(npz_fn, png_fn)

if __name__ == '__main__':
    main()
