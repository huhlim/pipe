#!/usr/bin/env python

import os
import sys
import path
import mdtraj
import numpy as np
from itertools import product

from libcommon import *

PARAM_DIST_SAME  = 0.2 # nm
PARAM_DIST_COORD = 0.3 # nm
PARAM_DIST_BSITE = 0.5 # nm
PARAM_DIST_ALIGN = 0.8 # nm

def get_ligand_info(job, wait_after_run, sleep=30):
    ligand_home = job.work_home.subdir("ligand", build=True)
    #
    with open("%s/STDRES"%DEFAULT_HOME) as fp:
        STDRES = []
        for line in fp:
            STDRES.append(line.strip())
    with open("%s/METALs"%DEFAULT_HOME) as fp:
        METALs = []
        for line in fp:
            METALs.append(line.strip())
    #
    str_fn_s = []
    while True:
        status = True
        #
        pdb_fn = ligand_home.glob("*.pdb")
        if len(pdb_fn) == 0:
            status = False
        ligand_pdb_fn = pdb_fn[0]
        #
        ligand_s = read_ligand_pdb(ligand_pdb_fn, STDRES, METALs)
        for ligName, is_metal in ligand_s:
            if is_metal: continue
            str_fn = ligand_home.fn("%s.str"%ligName)
            if not str_fn.status():
                status = False ; break
            str_fn_s.append(str_fn)
        #
        if wait_after_run:
            sys.stderr.write("waiting for ligand info... \n")
            time.sleep(sleep)
        else:
            break
    #
    set_bsite(job, ligand_pdb_fn)
    #
    if status: 
        job.ligand_pdb_fn = ligand_pdb_fn
        job.ligand_s = ligand_s
        job.str_fn_s = str_fn_s
        job.to_json()
    return status

def read_ligand_pdb(fn, STDRES, METALs):
    ligand_s = []
    with fn.open() as fp:
        for line in fp:
            if (not line.startswith("ATOM")) and (not line.startswith("HETATM")):
                continue
            resName = line[17:21].strip()
            if resName not in STDRES:
                ligand = (resName, resName in METALs)
                if ligand not in ligand_s:
                    ligand_s.append(ligand)
    return ligand_s

def get_aligned_residues(ref, pdb):
    dr = pdb[None,:] - ref[:,None]
    dmtx = np.sqrt(np.sum(dr**2, -1))
    aligned = np.argmin(dmtx, -1)
    status = (np.min(dmtx, -1) > PARAM_DIST_SAME)
    aligned[status] = -1
    return aligned

def get_bsite_geometry(dist_s, atomIndex, cbIndex, caIndex):
    d_sc = dist_s[atomIndex]
    d_cb = dist_s[cbIndex]
    d_ca = dist_s[caIndex]
    #
    geom_sc = []
    for i,j in zip(*np.where(d_sc < PARAM_DIST_COORD)):
        i_atm = atomIndex[i]
        d = d_sc[i,j]
        geom_sc.append((i_atm, j, d))
    if len(geom_sc) == 0:
        d = np.min(d_sc)
        i,j = np.where(d_sc == d)
        i_atm = atomIndex[i[0]]
        j = j[0]
        geom_sc.append((i_atm, j, d))
    #
    geom_ca = [(caIndex[0], np.argmin(d_ca), np.min(d_ca))]
    if cbIndex.shape[0] != 0:
        geom_cb = [(cbIndex[0], np.argmin(d_cb), np.min(d_cb))]
    else:
        geom_cb = None
    return geom_sc, geom_cb, geom_ca

def set_bsite(job, ligand_pdb_fn):
    ligand_pdb = mdtraj.load(ligand_pdb_fn.short())
    ligand_pdb = ligand_pdb.atom_slice(ligand_pdb.top.select("element != H"))
    #
    proteinIndex = ligand_pdb.top.select("protein")
    ligandIndex = ligand_pdb.top.select("not protein")
    #
    protein = ligand_pdb.atom_slice(proteinIndex)
    ligand = ligand_pdb.atom_slice(ligandIndex)
    ligandResidueNumber_s = np.unique([r.index for r in ligand.top.residues]).tolist()
    #
    atomPair_s = [atomPair for atomPair in product(proteinIndex, ligandIndex)]
    dist_s = mdtraj.compute_distances(ligand_pdb, atomPair_s)[0].reshape((len(proteinIndex), len(ligandIndex)))
    #
    bsiteAtom_s = proteinIndex[np.where(dist_s < PARAM_DIST_BSITE)[0]]
    bsiteResidue_s = np.unique([ligand_pdb.top.atom(i_atm).residue.index for i_atm in bsiteAtom_s])
    bsiteCalphaIndex = ligand_pdb.top.select("name CA")[bsiteResidue_s]
    bsiteCa_s = ligand_pdb.xyz[0,bsiteCalphaIndex]
    bsiteGeom_s = []
    for i_res in bsiteResidue_s:
        bsiteGeom_s.append(get_bsite_geometry(dist_s, \
                ligand_pdb.top.select("resi %d"%i_res), \
                ligand_pdb.top.select("resi %d and name CB"%i_res), \
                ligand_pdb.top.select("resi %d and name CA"%i_res)))
    #
    alignAtom_s = proteinIndex[np.where(dist_s < PARAM_DIST_ALIGN)[0]]
    alignResidue_s = np.unique([ligand_pdb.top.atom(i_atm).residue.index for i_atm in alignAtom_s])
    alignCalphaIndex = ligand_pdb.top.select("name CA")[alignResidue_s]
    alignCa_s = ligand_pdb.xyz[0,alignCalphaIndex]
    #
    top = mdtraj.load(job.top_fn.short())
    topCalphaIndex = top.top.select("name CA")
    topCa_s = top.xyz[0,topCalphaIndex]
    #
    bsiteAligned = get_aligned_residues(bsiteCa_s, topCa_s)
    bsiteAlignedCalphaIndex = topCalphaIndex[bsiteAligned[bsiteAligned != -1]]
    bsiteAlignedResidue_s = [top.top.atom(i_atm).residue for i_atm in bsiteAlignedCalphaIndex]
    #
    job.bsite_geometry = [] ; j = -1
    for i,align in enumerate(bsiteAligned):
        if align == -1: continue
        j += 1
        #
        bsiteResName = ligand_pdb.top.residue(bsiteResidue_s[i]).name
        bsiteAlignedResName = bsiteAlignedResidue_s[j].name
        #
        if bsiteAlignedResName == bsiteResName:  # use geom_sc
            geom = bsiteGeom_s[i][0]
        elif bsiteAlignedResName == 'GLY' or bsiteResName == 'GLY':  # use geom_ca
            geom = bsiteGeom_s[i][2]
        else:   # use either geom_ca or geom_cb depending on distance
            if bsiteGeom_s[i][1][0][2] < bsiteGeom_s[i][2][0][2]:
                geom = bsiteGeom_s[i][1]
            else:
                geom = bsiteGeom_s[i][2]
        #
        for i_prot, i_lig, d in geom:
            protResidue = bsiteAlignedResidue_s[j].index
            protAtomName = protein.top.atom(i_prot).name
            ligAtom = ligand.top.atom(i_lig)
            ligResidueIndex = ligandResidueNumber_s.index(ligAtom.residue.index)
            ligAtomName = ligAtom.residue.name
            job.bsite_geometry.append((protResidue, protAtomName, ligResidueIndex, ligAtomName, float(d)))
    #
    alignAligned = get_aligned_residues(alignCa_s, topCa_s)
    alignAlignedCalphaIndex = topCalphaIndex[alignAligned[alignAligned != -1]]
    alignAlignedResidue_s = [top.top.atom(i_atm).residue for i_atm in alignAlignedCalphaIndex]
    #
    job.bsite_align = [[], []] ; j = -1
    for i,align in enumerate(alignAligned):
        if align == -1: continue
        j += 1
        #
        alignResidueIndex = ligand_pdb.top.residue(alignResidue_s[i]).index
        alignAlignedResidueIndex = alignAlignedResidue_s[j].index
        job.bsite_align[0].append(alignResidueIndex)
        job.bsite_align[1].append(alignAlignedResidueIndex)

def add_ligand(out_fn, job, pdb_fn):
    ligand_pdb = mdtraj.load(job.ligand_pdb_fn.short())
    calphaIndex = ligand_pdb.top.select("name CA")
    ligand_align = calphaIndex[job.bsite_align[0]]
    #
    target_pdb = mdtraj.load(pdb_fn.short())
    calphaIndex = target_pdb.top.select("name CA")
    target_align = calphaIndex[job.bsite_align[1]]
    #
    ligand_pdb.superpose(target_pdb, atom_indices=ligand_align, ref_atom_indices=target_align)
    ligand = ligand_pdb.atom_slice(ligand_pdb.top.select("not protein"))
    for i,residue in enumerate(ligand.top.residues):
        residue.resSeq = i+1
    #
    top = target_pdb.top.join(ligand.top)
    xyz = np.concatenate([target_pdb.xyz[0], ligand.xyz[0]])[None,:]
    #
    pdb = mdtraj.Trajectory(xyz, top)
    pdb.save(out_fn.short())
    #
    het_s = {het[0][:3].strip(): het[0].strip() for het in job.ligand_s if len(het[0].strip()) > 3}
    if len(het_s) == 0:
        return
    #
    wrt = []
    with out_fn.open() as fp:
        for line in fp:
            if (not line.startswith("ATOM")) and (not line.startswith("HETATM")):
                wrt.append(line)
                continue
            #
            resName3 = line[17:21].strip()
            if resName3 not in het_s:
                wrt.append(line)
                continue
            #
            resName = het_s[resName3]
            line = '%s%4s%s'%(line[:17], resName, line[21:])
            wrt.append(line)
    with out_fn.open('wt') as fout:
        fout.writelines(wrt)

