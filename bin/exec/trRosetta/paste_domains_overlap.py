#!/usr/bin/env python

import os
import sys
import numpy as np
import argparse
import mdtraj

from libtrRosetta import system 

PARAM_BFACTOR_MAX = 99.99
PARAM_SIGNAL_RESIDUE = 5

def load_pdbs(fn_s):
    pdb_s = []
    for fn in fn_s:
        pdb_s.append(mdtraj.load(fn))
    return pdb_s

def read_bfactors(fn):
    bfac_s = {}
    with open(fn) as fp:
        for line in fp:
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                continue
            resNo = int(line[22:26])
            if resNo not in bfac_s:
                bfac_s[resNo] = []
            bfac = float(line[60:66])
            bfac_s[resNo].append(bfac)
    return bfac_s

def superpose_pdbs(init_pdb, refined_s):
    pdb_fn_s = [init_pdb] + refined_s
    pdb_s = load_pdbs(pdb_fn_s)
    resNo_s = []
    for pdb in pdb_s:
        resNo_s.append([residue.resSeq for residue in pdb.top.residues])
    resNo_refined = np.unique(np.concatenate(resNo_s[1:])).tolist()
    #
    c_exp = (resNo_s[0][-1] not in resNo_refined)
    if c_exp:
        resNo_s = [resNo_s[0]] + resNo_s[1:][::-1]
        pdb_s = [pdb_s[0]] + pdb_s[1:][::-1]
    #
    aligned = [True] + [False for _ in pdb_s[1:]]
    for i,ref in enumerate(pdb_s):
        for k,is_aligned in enumerate(aligned):
            if is_aligned: continue
            pdb = pdb_s[k]
            #
            overlap = np.intersect1d(resNo_s[i], resNo_s[k])
            ref_atoms = []
            for atom in ref.top.atoms:
                if atom.residue.resSeq in overlap and atom.name == 'CA':
                    ref_atoms.append(atom.index)
            pdb_atoms = []
            for atom in pdb.top.atoms:
                if atom.residue.resSeq in overlap and atom.name == 'CA':
                    pdb_atoms.append(atom.index)
            pdb.superpose(ref, atom_indices=pdb_atoms, ref_atom_indices=ref_atoms)
            aligned[k] = True
            if i == 0: break
    #
    if c_exp:
        pdb_s = [pdb_s[0]] + pdb_s[1:][::-1]

    return pdb_s[1:]

def paste_domains(init_pdb, refined_s, pdb_s):
    n_domain0 = len(pdb_s)
    pdb0 = mdtraj.load(init_pdb)
    residue_all = np.unique([r.resSeq for r in pdb0.top.residues])
    residue_all = np.sort(residue_all)
    #
    residue_s = []; xyz_s = []
    for pdb in pdb_s:
        residue_s.append([r.resSeq for r in pdb.top.residues])
        ca = pdb.atom_slice(pdb.top.select("name CA"))
        xyz_s.append(ca.xyz[0])
    #
    bfac_s = []
    for refined in refined_s:
        bfac_s.append(read_bfactors(refined))
    #
    domain = np.zeros_like(residue_all, dtype=int) - 1
    for i,residue in enumerate(residue_s):
        for resNo in residue:
            domain[residue_all==resNo] = i
    # residues which were not subjected to refine has -1.
    #
    if domain[0] == -1: # N-ter
        for i,resNo in enumerate(residue_all):
            if domain[i] == -1:
                resNo_nter = resNo
            else:
                break
        resNo_nter += PARAM_SIGNAL_RESIDUE
        #
        pdb = pdb0.atom_slice(pdb0.top.select("resSeq <= %d"%resNo_nter))
        pdb_s.append(pdb)
        #
        residue_s.append([r.resSeq for r in pdb.top.residues])
        ca = pdb.atom_slice(pdb.top.select("name CA"))
        xyz_s.append(ca.xyz[0])
        domain[residue_all <= resNo_nter] = len(residue_s)-1
        #
        bfac = {}
        for res in pdb.top.residues:
            bfac[res.resSeq] = np.ones(res.n_atoms) * PARAM_BFACTOR_MAX
        bfac_s.append(bfac)

    #
    if domain[-1] == -1: # C-ter
        for i,resNo in enumerate(residue_all[::-1]):
            if domain[-i-1] == -1:
                resNo_cter = resNo
            else:
                break
        resNo_cter -= PARAM_SIGNAL_RESIDUE
        #
        pdb = pdb0.atom_slice(pdb0.top.select("resSeq >= %d"%resNo_cter))
        pdb_s.append(pdb)
        #
        residue_s.append([r.resSeq for r in pdb.top.residues])
        ca = pdb.atom_slice(pdb.top.select("name CA"))
        xyz_s.append(ca.xyz[0])
        domain[residue_all >= resNo_cter] = len(residue_s)-1
        #
        bfac = {}
        for res in pdb.top.residues:
            bfac[res.resSeq] = np.ones(res.n_atoms) * PARAM_BFACTOR_MAX
        bfac_s.append(bfac)
    #
    n_domain = len(pdb_s)
    #
    for i in range(n_domain-1):
        for j in range(i+1, n_domain):
            linker_residue = np.intersect1d(residue_s[i], residue_s[j])
            if linker_residue.shape[0] == 0:
                continue
            #
            linker_s = []
            resPrev = None
            for res in linker_residue:
                if resPrev is None or res != resPrev+1:
                    linker = []
                    linker_s.append(linker)
                linker.append(res)
                resPrev = res
            #
            for linker in linker_s:
                xyz_i = xyz_s[i][np.isin(residue_s[i], linker)]
                xyz_j = xyz_s[j][np.isin(residue_s[j], linker)]
                dr_ij = xyz_j-xyz_i
                d_ij = np.sqrt(np.sum(dr_ij**2, axis=-1))
                center = linker[np.argmin(d_ij)]
                #
                if (linker[0]-1 in residue_s[i]) or (linker[-1]+1 in residue_s[j]):
                    # i: XXXXXXCX-----
                    # j: -----XCXXXXXX
                    I = i ; J = j
                else:
                    I = j ; J = i
                for resNo in linker:
                    if resNo < center:
                        domain[residue_all==resNo] = I
                    elif resNo > center:
                        domain[residue_all==resNo] = J
                    else:
                        domain[residue_all==resNo] = np.random.choice([I,J])
    #
    xyz = [] ; bfac = []
    top = mdtraj.Topology()
    chain = top.add_chain()
    for resNo,d in zip(residue_all, domain):
        residue = pdb_s[d].atom_slice(pdb_s[d].top.select("resSeq %d"%resNo))
        xyz.append(residue.xyz[0])
        bfac.append(bfac_s[d][resNo])
        #
        r = top.add_residue(residue.top.residue(0).name, chain, resSeq=resNo)
        for atom in residue.top.atoms:
            a = top.add_atom(atom.name, atom.element, r)

    xyz = np.concatenate(xyz)
    bfac = np.concatenate(bfac)
    #
    pdb = mdtraj.Trajectory(xyz[None,:], top)
    pdb.bfac = bfac
    return pdb

def main():
    arg = argparse.ArgumentParser(prog='average')
    arg.add_argument('--init', dest='init_pdb', required=True)
    arg.add_argument('--refined', dest='refined_pdbs', nargs='*', required=True)
    arg.add_argument('--output', dest='output_fn', required=True)
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    #
    pdb_s = superpose_pdbs(arg.init_pdb, arg.refined_pdbs)
    out = paste_domains(arg.init_pdb, arg.refined_pdbs, pdb_s)
    out.save(arg.output_fn, bfactors=out.bfac)

if __name__ == '__main__':
    main()
