#!/usr/bin/env python

import os
import sys
import mdtraj
import path
import numpy as np
import pickle
import argparse

from libanalysis import get_molecules, calc_center_of_mass, calc_average_boxsize, calc_distance


def get_random_distr(pbc, N=1000000):
    box = np.random.random(size=(N, 3)) * pbc - pbc / 2.0
    d = np.sqrt(np.sum(box**2, axis=1))
    return d


def get_histogram(d, dist_bin):
    h = np.histogram(d, bins=dist_bin)[0].astype(float)
    return h


def calc_distance_distr(cntr_i, cntr_j, pbc, dist_bin, is_same):
    pair_s = []
    for i, c_i in enumerate(cntr_i):
        for j, c_j in enumerate(cntr_j):
            if is_same and (j <= i):
                continue
            pair_s.append((i, j))
    if len(pair_s) == 0:
        return None
    #
    distr = np.zeros(dist_bin.shape[0] - 1, dtype=float)
    for (i, j) in pair_s:
        d = calc_distance(cntr_i[i], cntr_j[j], pbc)
        distr += get_histogram(d, dist_bin)
    return distr


def run(
    molecule_s, top, dcd_fn_s, T_STRIDE, T_SKIP, D_GRID, include_s=[], exclude_s=[], overwrite=False
):
    # get center of mass
    cntr_s = {}
    pbc = []
    for dcd_fn in dcd_fn_s:
        dcd_fn = path.Path(dcd_fn)
        cntr_fn = dcd_fn.dirname().fn(f"{dcd_fn.name()}.cntr.pkl")
        if cntr_fn.status() and (not overwrite):
            with cntr_fn.open("rb") as fp:
                X = pickle.load(fp)
            for molecule_name in X["cntr"]:
                if molecule_name in exclude_s:
                    continue
                if len(include_s) != 0 and molecule_name not in include_s:
                    continue
                if molecule_name not in cntr_s:
                    cntr_s[molecule_name] = []
                cntr_s[molecule_name].append(X["cntr"][molecule_name])
            pbc.append(X["pbc"])
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
        cntr = np.concatenate(cntr_s[molecule_name], axis=-2)[:, T_SKIP:]
        cntr_s[molecule_name] = cntr
    #
    pbc_max = pbc.max()
    pbc_avrg = pbc.mean()
    max_dist = np.ceil(pbc_max * (np.sqrt(3) / 2 / D_GRID)) * D_GRID
    dist_bin = np.arange(0, max_dist + 0.1 * D_GRID, D_GRID)
    #
    random_distr = get_random_distr(pbc_avrg)
    random_distr = get_histogram(random_distr, dist_bin)
    random_distr /= random_distr.sum()

    distr_s = {}
    distr_s["dist_bin"] = dist_bin * 10.0  # angstrom
    distr_s["random"] = random_distr
    #
    for i, mol_i in enumerate(cntr_s):
        for j, mol_j in enumerate(cntr_s):
            if j < i:
                continue
            #
            mol_pair = (mol_i, mol_j)
            d = calc_distance_distr(cntr_s[mol_i], cntr_s[mol_j], pbc, dist_bin, (i == j))
            if d is not None:
                distr_s[mol_pair] = d
    return distr_s


def main():
    arg = argparse.ArgumentParser(prog="calc_g_r")
    arg.add_argument(dest="out_fn")
    arg.add_argument("--top", dest="top_fn", required=True)
    arg.add_argument("--traj", dest="traj_s", required=True, nargs="*")
    arg.add_argument("--time_step", dest="time_step", default=0.1, type=float)
    arg.add_argument("--stride", dest="stride", default=1, type=int)
    arg.add_argument("--skip", dest="skip_first", default=0.0, type=float)
    arg.add_argument("--dgrid", dest="d_grid", default=0.5, type=float)
    arg.add_argument("--max_water", dest="max_water", default=-1, type=int)
    arg.add_argument("--include", dest="include_s", default=[], nargs="*")
    arg.add_argument("--exclude", dest="exclude_s", default=[], nargs="*")
    arg.add_argument("--overwrite", dest="overwrite", action="store_true", default=False)
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    #
    T_STRIDE = arg.stride
    T_SKIP = int(arg.skip_first // arg.time_step)
    D_GRID = arg.d_grid / 10.0  # Ang. -> nm
    #
    top = mdtraj.load(arg.top_fn)
    molecule_s = get_molecules(top, exclude_s=arg.exclude_s, max_water=arg.max_water)
    #
    distr_s = {}
    for traj_fn in arg.traj_s:
        try:
            with open(traj_fn) as fp:
                dcd_fn_s = [line.strip() for line in fp]
        except:
            dcd_fn_s = [traj_fn]
        #
        distr = run(
            molecule_s,
            top,
            dcd_fn_s,
            T_STRIDE,
            T_SKIP,
            D_GRID,
            include_s=arg.include_s,
            exclude_s=arg.exclude_s,
            overwrite=arg.overwrite,
        )
        for key in distr:
            if key not in distr_s:
                distr_s[key] = []
            distr_s[key].append(distr[key])
    #
    with open(arg.out_fn, "wb") as fout:
        pickle.dump(distr_s, fout)


if __name__ == "__main__":
    main()
