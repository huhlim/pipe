#!/usr/bin/env python

import os
import sys
import argparse
import tempfile
import numpy as np
from supPDB import supPDB
from libquat import sample_fibonacci_sphere, Quaternion

import mdtraj
from mdtraj.core.element import oxygen, hydrogen

WATER_NAME = ["HOH", "TIP", "WAT"]
WATER_OXYGEN = ["O", "OH2"]
WATER_GEOM = [0.09572, np.deg2rad(104.52)]
WATER_MAX_DISTANCE = 0.5
N_SAMPLE = 100
WATER_ORIENT = sample_fibonacci_sphere(n_sample=N_SAMPLE)


def select_water(pdblines, bfactor_cutoff=-1):
    wrt = []
    wat_s = []
    resNo_prev = None
    for line in pdblines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            resName = line[17:20].strip()
            if resName in WATER_NAME:
                resNo = line[21:27]
                if resNo != resNo_prev:
                    wat = [-1.0, []]
                    wat_s.append(wat)
                #
                atmName = line[12:16].strip()
                if atmName in WATER_OXYGEN:
                    try:
                        bfactor = float(line[60:68])
                    except:
                        bfactor = 0.0
                    wat[0] = bfactor
                    wat[1].append(line)
        elif line.startswith("END"):
            break
    #
    wat_s.sort(key=lambda x: [0])
    #
    for w in wat_s:
        if bfactor_cutoff < 0.0 or w[0] < bfactor_cutoff:
            wrt.extend(w[1])
    #
    water_pdb = tempfile.NamedTemporaryFile("wt", suffix=".pdb")
    water_pdb.writelines(wrt)
    #
    pdb = mdtraj.load(water_pdb.name)
    return pdb


def calc_min_d(xyz, r, return_index=False):
    dr = xyz[:, None] - r[None, :]
    d = np.sqrt(np.sum(dr**2, axis=-1))
    min_d = np.min(d)
    if return_index:
        return min_d, np.where(d == min_d)
    else:
        return min_d


def place_water(r_o, xyz, water_cutoff):
    r_h1 = r_o + WATER_ORIENT * WATER_GEOM[0]
    min_d, min_index = calc_min_d(xyz, r_h1, return_index=True)
    min_index = min_index[1][0]
    r_h1 = r_h1[min_index]
    v_h1 = WATER_ORIENT[min_index]
    #
    rv = np.random.randn(3)
    rv /= np.linalg.norm(rv)
    nv = np.cross(v_h1, rv)
    nv /= np.linalg.norm(nv)
    q0 = Quaternion.from_axis_and_angle(nv, WATER_GEOM[1])
    v_h2 = v_h1.dot(q0.rotate().T)
    r_h2 = []
    for angle in np.linspace(0.0, 2.0 * np.pi, N_SAMPLE):
        q = Quaternion.from_axis_and_angle(v_h1, angle)
        r_h2.append(r_o + v_h2.dot(q.rotate().T) * WATER_GEOM[0])
    r_h2 = np.array(r_h2)
    min_d, min_index = calc_min_d(xyz, r_h2, return_index=True)
    min_index = min_index[1][0]
    r_h2 = r_h2[min_index]
    #
    xyz = np.concatenate([xyz, np.array([r_o, r_h1, r_h2], dtype=np.float)])
    return xyz


def append_water(in_pdb, water, water_cutoff=0.2):
    pdb = mdtraj.load(in_pdb)
    top = pdb.topology
    xyz = pdb.xyz[0]
    #
    ia = 0
    chain = top.add_chain()
    resNo = 0
    water_index = 0
    segName = "WAT%1d" % water_index
    r_water = []
    for r, water in zip(water.xyz[0], water.topology.residues):
        min_d = calc_min_d(xyz, r)
        if min_d < water_cutoff:
            continue
        if min_d > WATER_MAX_DISTANCE:
            continue
        #
        if chain.n_residues >= 9999:
            chain = top.add_chain()
            resNo = 0
            water_index += 1
            segName = "WAT%1d" % water_index
        #
        resNo += 1
        residue = top.add_residue("TIP3", chain, resSeq=resNo, segment_id=segName)
        atom_o = top.add_atom("OH2", oxygen, residue)
        atom_h1 = top.add_atom("H1", hydrogen, residue)
        atom_h2 = top.add_atom("H2", hydrogen, residue)
        top.add_bond(atom_o, atom_h1)
        top.add_bond(atom_o, atom_h2)
        #
        xyz = place_water(r, xyz, water_cutoff)

    pdb = mdtraj.Trajectory(
        xyz[None, :],
        top,
        unitcell_lengths=pdb.unitcell_lengths,
        unitcell_angles=pdb.unitcell_angles,
    )
    return pdb


def main():
    arg = argparse.ArgumentParser(prog="convert_figure")
    arg.add_argument(dest="in_pdb")
    arg.add_argument(dest="ref_pdb")
    arg.add_argument(dest="out_pdb")
    arg.add_argument("-b", "--bfactor", dest="bfactor_cutoff", default=-1.0, type=float)
    arg.add_argument("-c", "--cutoff", dest="water_cutoff", default=1.0, type=float)
    #
    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args()
    #
    sup = supPDB(arg.ref_pdb)
    sup.align(arg.in_pdb, method="TMalign")
    water = select_water(sup.write(), bfactor_cutoff=arg.bfactor_cutoff)
    pdb = append_water(arg.in_pdb, water, water_cutoff=0.1 * arg.water_cutoff)
    pdb.save(arg.out_pdb)


if __name__ == "__main__":
    main()
