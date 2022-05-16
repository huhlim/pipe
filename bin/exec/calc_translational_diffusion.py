#!/usr/bin/env python

import os
import sys
import mdtraj
import path
import numpy as np
import pickle
import argparse

from libanalysis import get_molecules, calc_center_of_mass, calc_average_boxsize, unwrap_PBC

def run(molecule_s, top, dcd_fn_s, T_STEP, T_STRIDE, T_SKIP, T_MAX, T_LAG_STEP, \
        include_s=[], exclude_s=[], overwrite=False):
    cntr_s = {}
    pbc = []
    for dcd_fn in dcd_fn_s:
        dcd_fn = path.Path(dcd_fn)
        cntr_fn = dcd_fn.dirname().fn(f'{dcd_fn.name()}.cntr.pkl')
        if cntr_fn.status() and (not overwrite):
            with cntr_fn.open("rb") as fp:
                X = pickle.load(fp)
            for molecule_name in X['cntr']:
                if molecule_name not in cntr_s:
                    cntr_s[molecule_name] = []
                cntr_s[molecule_name].append(X['cntr'][molecule_name])
            pbc.append(X['pbc'])
        else:
            traj = mdtraj.load(dcd_fn.path(), top=top, stride=T_STRIDE)
            cntr = calc_center_of_mass(molecule_s, traj, include_s=include_s, exclude_s=exclude_s)
            for molecule_name in cntr:
                if molecule_name not in cntr_s:
                    cntr_s[molecule_name] = []
                cntr_s[molecule_name].append(cntr[molecule_name])
            pbc.append(traj.unitcell_lengths)
            #
            with cntr_fn.open("wb") as fout:
                pickle.dump({"cntr": cntr, "pbc": traj.unitcell_lengths}, fout)
        #
    pbc = np.concatenate(pbc, axis=0)[T_SKIP:]
    for molecule_name in cntr_s:
        cntr = np.concatenate(cntr_s[molecule_name], axis=-2)[:,T_SKIP:]
        cntr_s[molecule_name] = unwrap_PBC(cntr, pbc)
    #
    T_LAG = np.array([int(t/T_STEP) for t in np.arange(T_LAG_STEP, T_MAX+T_LAG_STEP, T_LAG_STEP)], dtype=int)
    T_LAG = np.unique(T_LAG)
    t_lag = T_LAG.astype(float) * T_STEP
    #
    msd_s = {}
    for molecule_name in cntr_s:
        cntr = cntr_s[molecule_name]
        msd = np.zeros((cntr.shape[0], len(T_LAG)), dtype=float)
        #
        for t,tau in enumerate(T_LAG):
            dr = cntr[:,:-tau] - cntr[:,tau:]
            dsq = np.sum(dr**2, axis=-1)
            msd[:,t] = np.mean(dsq, axis=-1)
        msd_s[molecule_name] = (t_lag, msd)
    #
    diffusion_s = {}
    for molecule_name in msd_s:
        diffusion = np.polyfit(t_lag, msd_s[molecule_name][1].T, 1)[0] / 6.
        diffusion_s[molecule_name] = diffusion

    return msd_s, diffusion_s

def main():
    arg = argparse.ArgumentParser(prog='calc_translational_diffusion')
    arg.add_argument(dest='out_fn')
    arg.add_argument('--top', dest="top_fn", required=True)
    arg.add_argument('--traj', dest='traj_s', required=True, nargs='*')
    arg.add_argument('--time_step', dest='time_step', required=True, type=float)
    arg.add_argument('--stride', dest='stride', default=1, type=int)
    arg.add_argument('--skip', dest='skip_first', default=0., type=float)
    arg.add_argument('--lag_time', dest='lag_time', default=1., type=float)
    arg.add_argument('--max_time', dest='max_time', default=100., type=float)
    arg.add_argument('--max_water', dest='max_water', default=-1, type=int)
    arg.add_argument('--include', dest='include_s', default=[], nargs='?')
    arg.add_argument('--exclude', dest='exclude_s', default=[], nargs='?')
    arg.add_argument('--overwrite', dest='overwrite', action='store_true', default=False)
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    #
    T_STEP = arg.time_step  # in ns
    T_STRIDE = arg.stride
    T_SKIP = int(arg.skip_first // T_STEP)
    T_MAX = arg.max_time
    T_LAG_STEP = arg.lag_time
    #
    top = mdtraj.load(arg.top_fn)
    molecule_s = get_molecules(top, exclude_s=arg.exclude_s, max_water=arg.max_water)
    #
    msd_s = {}
    diffusion_s = {}
    for traj_fn in arg.traj_s:
        try:
            with open(traj_fn) as fp:
                dcd_fn_s = [line.strip() for line in fp]
        except:
            dcd_fn_s = [traj_fn]
        #
        msd, diffusion = run(molecule_s, top, dcd_fn_s, T_STEP, T_STRIDE, T_SKIP, T_MAX, T_LAG_STEP, \
                include_s=arg.include_s, exclude_s=arg.exclude_s, overwrite=arg.overwrite)
        #
        for molecule_name in diffusion:
            if molecule_name not in msd_s:
                msd_s[molecule_name] = []
                diffusion_s[molecule_name] = []
            #
            msd_s[molecule_name].append(msd[molecule_name])
            diffusion_s[molecule_name].append(diffusion[molecule_name])
    #
    with open(arg.out_fn, 'wb') as fout:
        pickle.dump({"msd": msd_s, "diffusion": diffusion_s}, fout)


if __name__ == '__main__':
    main()

