#!/usr/bin/env python

import os
import sys
import mdtraj
import argparse
import numpy as np

def main():
    arg = argparse.ArgumentParser(prog='prod_meta')
    arg.add_argument('--index', dest='index_fn', required=True)
    arg.add_argument('--top', dest='top_fn', required=True)
    arg.add_argument('--dcd', dest='dcd_fn_s', nargs='*', default=[])
    arg.add_argument('--out', dest='out_fn', required=True)
    #
    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args()
    if len(arg.dcd_fn_s) == 0:
        return
    #
    top = mdtraj.load(arg.top_fn)
    atomIndex = np.load(arg.index_fn)['select']
    #
    traj_s = []
    for dcd_fn in arg.dcd_fn_s:
        traj = mdtraj.load(dcd_fn, top=top, atom_indices=atomIndex)
        traj_s.append(traj)
    traj_s = mdtraj.join(traj_s, check_topology=False)
    traj_s.save(arg.out_fn)

if __name__ == '__main__':
    main()
