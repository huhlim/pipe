#!/usr/bin/env python

import os
import sys
import copy
import numpy as np

from libtrRosetta import *
from get_domain import get_soft_domain_boundary

sys.path.insert(0, '../../')
import path

APP_HOME = "/home/huhlim/apps/trRosetta"
PARAM_CUTOFF = 0.15
FEATUREs = ['dist', 'omega', 'theta', 'phi']

def run(msa_fn):
    work_home = msa_fn.dirname()
    work_home.chdir()
    #
    title = msa_fn.name()
    #
    out_fn = work_home.fn("%s.npz"%title)
    if out_fn.status():
        return out_fn
    #
    cmd = []
    cmd.append("python")
    cmd.append("%s/network/predict.py"%APP_HOME)
    cmd.extend(['-m', "%s/model2019_07"%APP_HOME])
    cmd.append(msa_fn.short())
    cmd.append(out_fn.short())
    #
    log_fn = work_home.fn("trRosetta.log")
    with log_fn.open("wt") as fout:
        system(cmd, outfile=fout, redirect=True)
    #
    return out_fn

def get_acc(dist):
    conf = dist[:,:,1:13].sum(axis=-1)
    l_seq = conf.shape[0]
    mlc_index = np.triu_indices(l_seq, k=12)
    mlc = conf[mlc_index]
    mlc.sort()
    mlc = mlc[::-1][:l_seq]
    return mlc.mean()

def merge(npz0, npz1, domain):
    if npz0 is None:
        npz1_acc = get_acc(npz1['dist'])
        sys.stdout.write("MERGE: %-12s  %6.4f %6.4f\n"%(domain.name, npz1_acc, npz1_acc))
        return npz1
    #
    mask = np.ix_(domain.mask, domain.mask)
    #
    npz0_local = npz0['dist'][mask]
    npz1_local = npz1['dist']
    npz0_local_acc = get_acc(npz0_local)
    npz1_local_acc = get_acc(npz1_local)
    if npz0_local_acc > npz1_local_acc:
        return npz0
    #
    n_domain_res = npz0_local.shape[0]
    w_1d = get_soft_domain_boundary(domain.mask)[domain.mask]
    w_2d = np.tile(w_1d, [n_domain_res, 1])
    w_2d = np.minimum(w_2d, w_2d.T)
    w_3d = w_2d[:,:,None]
    #
    npz2 = copy.deepcopy(npz0)
    for feature in FEATUREs:
        npz2[feature][mask] = w_3d * npz1[feature] + (1.-w_3d) * npz0[feature][mask]
    #
    npz0_global_acc = get_acc(npz0['dist'])
    npz2_global_acc = get_acc(npz2['dist'])
    if npz0_global_acc > npz2_global_acc:
        return npz0
    else:
        sys.stdout.write("MERGE: %-12s  %6.4f %6.4f -> %6.4f %6.4f\n"%(domain.name, \
                npz0_global_acc, npz0_local_acc, \
                npz2_global_acc, npz1_local_acc))
        return npz2

def get_contact_trRosetta(npz_fn):
    npz = np.load(npz_fn.short())
    prob = npz['dist'][:,:,1:13].sum(axis=-1)
    prob[prob < PARAM_CUTOFF] = 0.0
    return prob

def read_trRosetta(npz_fn):
    npz = np.load(npz_fn.short())
    out = {}
    for feature in FEATUREs:
        out[feature] = npz[feature]
    return out

def main():
    in_fn = path.Path(sys.argv[1])
    run(in_fn)

if __name__ == '__main__':
    main()
