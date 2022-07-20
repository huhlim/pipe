#!/usr/bin/env python

import os
import sys
import copy
import numpy as np
import mdtraj
import argparse

from libquat import Quaternion as Quat
from string import ascii_uppercase as CHAINs

N_AVOGADRO = 6.0221420e23
PARAM_CLASH_CUTOFF = 0.35  # [nm]
PARAM_MAX_TRY = 100


def calc_distsance_pbc_1(X, Y, BOX_WIDTH):
    BOX_HALF = BOX_WIDTH * 0.5
    dr = X - Y[:, np.newaxis]
    dr[dr > BOX_HALF] = BOX_WIDTH - dr[dr > BOX_HALF]
    dr[dr < -BOX_HALF] = BOX_WIDTH + dr[dr < -BOX_HALF]
    #
    d = np.sqrt(np.sum(dr**2, axis=-1))
    return d


def calc_distsance_pbc_2(X, Y, BOX_WIDTH):
    BOX_HALF = BOX_WIDTH * 0.5
    #
    d = np.zeros((Y.shape[0], X.shape[0]), dtype=np.float16)
    for i, y in enumerate(Y):
        dr = X - y
        dr[dr > BOX_HALF] = BOX_WIDTH - dr[dr > BOX_HALF]
        dr[dr < -BOX_HALF] = BOX_WIDTH + dr[dr < -BOX_HALF]
        d[i] = np.sqrt(np.sum(dr**2, axis=-1))
    return d


def check_clash(X, Y, BOX_WIDTH, cutoff):
    if len(Y) == 0:
        return False
    #
    rX = np.array(X).astype(np.float16)
    nX = rX.shape[0]
    rY = np.concatenate(Y).astype(np.float16)
    nY = rY.shape[0]
    #
    if nX * nY < 1000000000 or False:
        min_d = np.min(calc_distsance_pbc_1(rX, rY, BOX_WIDTH))
    else:
        min_d = np.min(calc_distsance_pbc_2(rX, rY, BOX_WIDTH))
    return min_d < cutoff


def place(molecule_s, n_place, box_width):
    xyz_s = []
    xyz_hv_s = []
    for molecule in molecule_s:
        xyz = molecule.xyz[0]
        xyz_hv = molecule.xyz[0, molecule.top.select("element != H")]
        cntr = mdtraj.compute_center_of_mass(molecule)[0]
        #
        xyz -= cntr
        xyz_s.append(xyz)
        xyz_hv -= cntr
        xyz_hv_s.append(xyz_hv)
    #
    n_mol_type = len(molecule_s)
    #
    def place_runner(clash_cutoff):
        system = []
        system_hv = []
        for i in range(n_mol_type):
            for k in range(n_place[i]):
                status = False
                for _ in range(PARAM_MAX_TRY):
                    tr = np.random.random(size=3) * box_width
                    q = np.random.random(size=4)
                    q /= np.sqrt(np.sum(q**2))
                    q = Quat(q)
                    #
                    r_hv = np.matmul(xyz_hv_s[i], q.rotate()) + tr
                    #
                    if not check_clash(r_hv, system_hv, box_width, clash_cutoff):
                        status = True
                        system_hv.append(r_hv)
                        r = np.matmul(xyz_s[i], q.rotate()) + tr
                        system.append(r)
                        break
                #
                if not status:
                    return False, None
        #
        return True, system

    #
    N_FAIL = 0
    while True:
        clash_cutoff = PARAM_CLASH_CUTOFF - 0.01 * int(N_FAIL / 100.0)
        status, system = place_runner(clash_cutoff)
        if status:
            return system
        else:
            N_FAIL += 1


def convert_to_pdb(molecule_s, n_place, xyz, box_width):
    xyz = np.concatenate([r.reshape((1, -1, 3)) for r in xyz], axis=1)
    #
    top = mdtraj.Topology()
    #
    for i, n_mol in enumerate(n_place):
        molecule = molecule_s[i]
        for k in range(n_place[i]):
            top = top.join(molecule.top)
    for i, chain in enumerate(top.chains):
        for residue in chain.residues:
            residue.segment_id = "P%03d" % i

    out = mdtraj.Trajectory(
        xyz, top, unitcell_lengths=box_width * np.ones(3), unitcell_angles=90.0 * np.ones(3)
    )
    return out


def run(in_fn_s, n_place, box_width, keep_atom_name):
    molecule_s = [mdtraj.load(fn, standard_names=(not keep_atom_name)) for fn in in_fn_s]
    xyz = place(molecule_s, n_place, box_width)
    out = convert_to_pdb(molecule_s, n_place, xyz, box_width)
    return out


def get_mass(fn):
    pdb = mdtraj.load(fn)
    mass = 0.0
    for atom in pdb.top.atoms:
        mass += atom.element.mass
    return mass


def calc_system_size(arg):
    n_molecules = len(arg.in_fn_s)
    n_place = np.array(arg.n_place, dtype=int)
    if len(arg.conc) == n_molecules:
        conc = np.array(arg.conc, dtype=float)
    else:
        conc = np.zeros(n_molecules, dtype=float)
    conc_total = arg.conc_total
    mass = np.array([get_mass(fn) for fn in arg.in_fn_s], dtype=float)
    if arg.unit == "g/L":
        conc /= mass * 0.001
    box_width = arg.box_width * 0.1  # nm
    box_volume = box_width**3  # nm^3 = 10^-24 L
    #
    status = True
    if np.all(n_place > 0):
        if box_width > 0.0:
            pass
        elif np.any(conc > 0.0):
            n = n_place[conc > 0.0]
            c = conc[conc > 0.0]
            v = n / N_AVOGADRO / (c * 1e-3)  # Liter
            box_width = (v * 1e-3) ** (1.0 / 3.0) / 1e-9  # nm
            box_width = box_width.mean()
        elif conc_total > 0.0:
            if arg.unit == "mM":
                n = n_place.astype(float).sum()
                v = n / N_AVOGADRO / (conc_total * 1e-3)  # Liter
            else:
                m = (n_place * mass).sum() / N_AVOGADRO  # gram
                v = m / conc_total  # Liter
            box_width = (v * 1e-3) ** (1.0 / 3.0) / 1e-9  # nm
        else:
            status = False
        box_volume = box_width**3  # nm^3 = 10^-24 L
        if status:
            conc = n_place.astype(float) / N_AVOGADRO / (box_volume * 1e-24) * 1e3
    elif np.any(n_place > 0):
        if np.any(conc[n_place > 0] > 0.0):
            subset = (n_place > 0) & (conc > 0.0)
            n = n_place[subset]
            c = conc[subset]
            v = n / N_AVOGADRO / (c * 1e-3)  # Liter
            box_width = (v * 1e-3) ** (1.0 / 3.0) / 1e-9  # nm
            box_width = box_width.mean()
            box_volume = box_width**3  # nm^3 = 10^-24 L

            subset = (n_place > 0) & (conc <= 0.0)
            conc[subset] = n_place[subset].astype(float) / N_AVOGADRO / box_volume * 1e3
            if np.all(conc > 0.0):
                n = (conc * 1e-3) * (box_volume * 1e-24) * N_AVOGADRO
                n_place[n_place <= 0] = np.round(n[n_place <= 0])
            elif np.where(conc <= 0.0)[0].shape[0] == 1 and conc_total > 0.0:
                if arg.unit == "mM":
                    c = conc[conc > 0.0].sum()
                    conc[conc <= 0.0] = conc_total - c
                else:
                    c = (conc[conc > 0.0] * 1e-3) * mass
                    c = conc_total - c.sum()
                    conc[conc <= 0.0] = c / mass[conc <= 0.0] * 1e3
                n_place = (conc * 1e-3) * (box_volume * 1e-24) * N_AVOGADRO
            else:
                status = False
        else:
            status = False
    else:
        if box_width > 0.0:
            if np.all(conc > 0.0):
                n_place = (conc * 1e-3) * (box_volume * 1e-24) * N_AVOGADRO
            elif np.where(conc <= 0.0)[0].shape[0] == 1 and conc_total > 0.0:
                if arg.unit == "mM":
                    c = conc[conc > 0.0].sum()
                    conc[conc <= 0.0] = conc_total - c
                else:
                    c = (conc[conc > 0.0] * 1e-3) * mass
                    c = conc_total - c.sum()
                    conc[conc <= 0.0] = c / mass[conc <= 0.0] * 1e3
                n_place = (conc * 1e-3) * (box_volume * 1e-24) * N_AVOGADRO
            else:
                status = False
        else:
            status = False

    if not status:
        sys.exit("ERROR: cannot determine the number of molecules.")
    n_place = n_place.astype(int)
    box_volume = box_width**3  # nm^3 = 10^-24 L
    #
    sys.stdout.write("# INPUT arguments\n")
    sys.stdout.write("BOX_WIDTH   %9.4f [nm]\n" % (arg.box_width * 0.1))
    sys.stdout.write("TOTAL_CONC  %9.4f [%s]\n" % (arg.conc_total, arg.unit))
    for i, in_fn in enumerate(arg.in_fn_s):
        wrt = []
        wrt.append("%5d" % arg.n_place[i])
        wrt.append("%9.4f [%s]" % (arg.conc[i], arg.unit))
        wrt.append("    %s" % in_fn)
        sys.stdout.write(" ".join(wrt) + "\n")
    sys.stdout.write("#\n")
    #
    sys.stdout.write("# FINAL compositions\n")
    sys.stdout.write("BOX_WIDTH   %9.4f [nm]\n" % box_width)
    sys.stdout.write(
        "TOTAL_CONC  %9.4f [mM]    %9.4f [g/L]\n" % (conc.sum(), ((conc * 1e-3) * mass).sum())
    )
    sys.stdout.write(
        "SYSTEM_SIZE   %7d [atoms]\n" % ((box_width * 1e-9) ** 3 / 1e-6 / 18.01528 * N_AVOGADRO * 3)
    )
    for i, in_fn in enumerate(arg.in_fn_s):
        wrt = []
        wrt.append("%5d" % n_place[i])
        wrt.append("%9.4f [mM]" % (conc[i]))
        wrt.append("%9.4f [g/L]" % (conc[i] * mass[i] * 1e-3))
        wrt.append("%9.2f [g/mol]" % mass[i])
        wrt.append("    %s" % in_fn)
        sys.stdout.write(" ".join(wrt) + "\n")
    return n_place, box_width


def main():
    arg = argparse.ArgumentParser(prog="build_system")
    arg.add_argument("out_fn")
    arg.add_argument("-i", "--input", dest="in_fn_s", required=True, nargs="*")
    arg.add_argument("-n", dest="n_place", nargs="*", type=int, default=[])
    arg.add_argument("-c", dest="conc", nargs="*", type=float, default=[])
    arg.add_argument("--conc", dest="conc_total", type=float, default=0.0)
    arg.add_argument("--unit", dest="unit", type=str, choices=["mM", "g/L"], default="mM")
    arg.add_argument("-b", "--box", dest="box_width", type=float, default=0.0)
    arg.add_argument("--check", dest="check_system_size", action="store_true", default=False)
    arg.add_argument("--keep_atom_name", dest="keep_atom_name", default=False, action="store_true")
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    if len(arg.n_place) > 0 and len(arg.n_place) != len(arg.in_fn_s):
        sys.exit("ERROR: number of molecules does not match.")
    if len(arg.conc) > 0 and len(arg.conc) != len(arg.in_fn_s):
        sys.exit("ERROR: number of conc does not match.")
    if len(arg.n_place) == 0:
        arg.n_place = [0 for _ in arg.in_fn_s]
    if len(arg.conc) == 0:
        arg.conc = [0.0 for _ in arg.in_fn_s]
    arg.in_fn_s = [os.path.abspath(fn) for fn in arg.in_fn_s]
    #
    n_place, box_width = calc_system_size(arg)
    if arg.check_system_size:
        return
    pdb = run(arg.in_fn_s, n_place, box_width, arg.keep_atom_name)
    pdb.save(arg.out_fn)


if __name__ == "__main__":
    main()
