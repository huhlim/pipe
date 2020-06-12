#!/usr/bin/env python

import os
import sys
import numpy as np
import mdtraj

from pyrosetta import *
from pyrosetta.rosetta.core import scoring
from pyrosetta.rosetta.protocols.minimization_packing import MinMover

def get_structure_info(pdb_fn):
    pdb = mdtraj.load(pdb_fn)
    calpha = pdb.atom_slice(pdb.top.select("name CA"))
    resNo_s = [r.resSeq for r in calpha.top.residues]
    #
    Rs = calpha.xyz[0]
    dR = Rs[None,:] - Rs[:,None]
    dij = np.sqrt(np.sum(dR**2, -1)) * 10.0 # angstrom
    l_seq = Rs.shape[0]
    valid = np.ones(l_seq, dtype=bool)
    valid[:3] = False
    valid[-3:] = False
    #
    i,j = np.where(dij < 20.0)
    dij = dij[i,j]
    cst_s = []
    for ii,jj,dd in zip(i,j,dij):
        if not (valid[ii] and valid[jj]):
            continue
        i_res = resNo_s[ii]
        j_res = resNo_s[jj]
        if j_res - i_res < 2:
            continue
        id_i = rosetta.core.id.AtomID(2, i_res)
        id_j = rosetta.core.id.AtomID(2, j_res)
        #
        harmonic = rosetta.core.scoring.func.HarmonicFunc(dd, 0.5)
        cst = rosetta.core.scoring.constraints.AtomPairConstraint(id_i, id_j, harmonic)
        cst_s.append(cst)
    #
    phi_s = mdtraj.compute_phi(pdb)[1][0]
    psi_s = mdtraj.compute_psi(pdb)[1][0]
    omg_s = mdtraj.compute_omega(pdb)[1][0]
    dihed_s = {}
    for i,resNo in enumerate(resNo_s):
        dihed_s[resNo] = np.concatenate([(2.*np.random.random(2)-1.)*np.pi, [np.pi]])
    for i,resNo in enumerate(resNo_s[1:]):
        if valid[resNo_s.index(resNo)]:
            dihed_s[resNo][0] = phi_s[i]
    for i,resNo in enumerate(resNo_s[:-1]):
        if valid[resNo_s.index(resNo)]:
            dihed_s[resNo][1] = psi_s[i]
            dihed_s[resNo][2] = omg_s[i]

    for i,resNo in enumerate(resNo_s):
        dihed_s[resNo] = np.degrees(dihed_s[resNo])
    return dihed_s, cst_s

def set_dihed(pose, dihed_s):
    for resNo in range(1, pose.total_residue()+1):
        if resNo in dihed_s:
            pose.set_phi(resNo, dihed_s[resNo][0])
            pose.set_psi(resNo, dihed_s[resNo][1])
            pose.set_omega(resNo, dihed_s[resNo][2])
        else:
            pose.set_phi(resNo, np.random.random()*360.0-180.0)
            pose.set_psi(resNo, np.random.random()*360.0-180.0)
            pose.set_omega(resNo, 180.)
    return pose

def set_scoring_function():
    sf_s = []
    #
    # sf_vdw
    sf = ScoreFunction()
    sf.set_weight(scoring.rama, 1.0)
    sf.set_weight(scoring.vdw, 1.0)
    sf_s.append(sf)
    #
    # sf_vdw2
    sf = ScoreFunction()
    sf.set_weight(scoring.rama, 1.0)
    sf.set_weight(scoring.vdw, 1.0)
    sf.set_weight(scoring.atom_pair_constraint, 5.0)
    sf_s.append(sf)
    #
    # sf
    sf = ScoreFunction()
    sf.set_weight(scoring.cen_hb, 5.0)
    sf.set_weight(scoring.rama, 1.0)
    sf.set_weight(scoring.omega, 0.5)
    sf.set_weight(scoring.vdw, 1.0)
    sf.set_weight(scoring.atom_pair_constraint, 5.0)
    sf.set_weight(scoring.dihedral_constraint, 4.0)
    sf.set_weight(scoring.angle_constraint, 4.0)
    sf_s.append(sf)
    #
    # sf1
    sf = ScoreFunction()
    sf.set_weight(scoring.cen_hb, 5.0)
    sf.set_weight(scoring.rama, 1.0)
    sf.set_weight(scoring.omega, 0.5)
    sf.set_weight(scoring.vdw, 3.0)
    sf.set_weight(scoring.atom_pair_constraint, 5.0)
    sf.set_weight(scoring.dihedral_constraint, 1.0)
    sf.set_weight(scoring.angle_constraint, 1.0)
    sf_s.append(sf)
    #
    # sf_vdw_cart
    sf = ScoreFunction()
    sf.set_weight(scoring.hbond_sr_bb, 3.0)
    sf.set_weight(scoring.hbond_lr_bb, 3.0)
    sf.set_weight(scoring.rama, 1.0)
    sf.set_weight(scoring.omega, 0.5)
    sf.set_weight(scoring.vdw, 0.5)
    sf.set_weight(scoring.cart_bonded, 0.1)
    sf.set_weight(scoring.atom_pair_constraint, 5.0)
    sf.set_weight(scoring.dihedral_constraint, 4.0)
    sf.set_weight(scoring.angle_constraint, 4.0)
    sf_s.append(sf)
    #
    sf = create_score_function('ref2015')
    sf.set_weight(scoring.atom_pair_constraint, 5.)
    sf.set_weight(scoring.dihedral_constraint, 1.)
    sf.set_weight(scoring.angle_constraint, 1.)
    sf_s.append(sf)

    return sf_s

def set_mover(sf_s):
    mmap = MoveMap()
    mmap.set_bb(True)
    mmap.set_chi(False)
    mmap.set_jump(True)
    #
    mover_s = []
    for sf in sf_s[:-1]:
        mover = MinMover(mmap, sf, 'lbfgs_armijo_nonmonotone', 0.0001, True)
        mover.max_iter(500)
        mover_s.append(mover)
    mover_s[-1].cartesian(True)
    #
    mmap = MoveMap()
    mmap.set_bb(True)
    mmap.set_chi(True)
    mmap.set_jump(True)

    relax = rosetta.protocols.relax.FastRelax()
    relax.set_scorefxn(sf_s[-1])
    relax.max_iter(100)
    relax.dualspace(True)
    relax.set_movemap(mmap)
    mover_s.append(relax)
    return mover_s

def remove_clash(sf, mover, pose):
    for _ in range(5):
        if float(sf(pose)) < 10.:
            break
        mover.apply(pose)

def run_pyRosetta(out_fn, sequence, pdb_fn_s):
    init('-mute all')
    #
    sf_s = set_scoring_function()
    mover_s = set_mover(sf_s)
    cset = rosetta.core.scoring.constraints.ConstraintSet()
    #
    dihed_s = {}
    for pdb_fn in pdb_fn_s:
        struct_s = get_structure_info(pdb_fn)
        dihed_s.update(struct_s[0])
        for cst in struct_s[1]:
            cset.add_constraint(cst)
    #
    pose = pose_from_sequence(sequence, 'centroid')
    set_dihed(pose, dihed_s)
    #
    cst = rosetta.protocols.constraint_movers.ConstraintSetMover()
    cst.constraint_set(cset)
    cst.add_constraints(True)
    cst.apply(pose)
    #
    remove_clash(sf_s[0], mover_s[1], pose)
    #
    mover_s[2].apply(pose)
    mover_s[4].apply(pose)
    remove_clash(sf_s[0], mover_s[3], pose)
    #
    switch = SwitchResidueTypeSetMover("fa_standard")
    switch.apply(pose)
    mover_s[-1].apply(pose)
    #
    pose.dump_pdb(out_fn)

def main():
    if len(sys.argv) == 1:
        sys.exit("usage: %s [OUT] [FA] [PDBs]"%__file__)
    #
    out_fn = sys.argv[1]
    in_fa = sys.argv[2]
    pdb_fn_s = sys.argv[3:]
    #
    with open(in_fa) as fp:
        sequence = []
        for line in fp:
            if not line.startswith(">"):
                sequence.append(line.strip())
        sequence = ''.join(sequence)
    #
    run_pyRosetta(out_fn, sequence, pdb_fn_s)

if __name__ == '__main__':
    main()
