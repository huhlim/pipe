#!/usr/bin/env python

import os
import sys
import path
import numpy as np
import mdtraj
from itertools import product
from run_tbm import METHODs

np.set_printoptions(suppress=True, linewidth=1000, edgeitems=100)

FEATUREs = ['dist', 'omega', 'theta', 'phi']

def get_distance(residue1, residue2, pdb_s):
    pair_s = product(residue1, residue2)
    dij = mdtraj.compute_distances(pdb_s, pair_s)
    dij = dij.reshape(-1, len(residue1), len(residue2))
    return dij*10.0

def get_angle(residue1, residue2, pdb_s):
    pair_s = product(residue1, residue2)
    pair_s = [np.concatenate(pair) for pair in pair_s]
    angle = mdtraj.compute_angles(pdb_s, pair_s)
    angle = angle.reshape(-1, len(residue1), len(residue2))
    return angle

def get_dihedral(residue1, residue2, pdb_s):
    pair_s = product(residue1, residue2)
    pair_s = [np.concatenate(pair) for pair in pair_s]
    angle = mdtraj.compute_dihedrals(pdb_s, pair_s)
    angle = angle.reshape(-1, len(residue1), len(residue2))
    return angle

def get_features(pdb_s):
    top = pdb_s.topology
    is_gly = np.array([(residue.name == 'GLY') for residue in top.residues], dtype=bool)
    #
    atomIndex = np.zeros((4, top.n_residues), dtype=int) # Ca, Cb, N, O
    atomIndex[0] = top.select("name CA")
    atomIndex[1,~is_gly] = top.select("name CB")
    atomIndex[1, is_gly] = top.select("name CA")[is_gly]
    atomIndex[2] = top.select("name N")
    atomIndex[3] = top.select("name O")
    #
    feature_s = {}
    #
    feature_s['dist']   = get_distance(atomIndex[1], atomIndex[1], pdb_s)
    feature_s['phi']    = get_angle(atomIndex[(0,1),].T, atomIndex[(1,),].T, pdb_s)
    feature_s['omega']  = get_dihedral(atomIndex[(0,1),].T, atomIndex[(1,0),].T, pdb_s)
    feature_s['theta']  = get_dihedral(atomIndex[(2,0,1),].T, atomIndex[(1,),].T, pdb_s)
    #
    return is_gly, feature_s

def get_hist(n_residue, feature_s):
    DIST_BINS = np.array(np.arange(2.0, 20.0, 0.5).tolist() + [20.0, 9999.9])
    ANGLE_BINS = np.deg2rad(np.arange(0.0, 180.0, 15.0))
    DIHEDRAL_BINS = np.deg2rad(np.arange(0.0, 360.0, 15.0))
    BINS = {'dist': DIST_BINS, 'phi': ANGLE_BINS, 'omega': DIHEDRAL_BINS, 'theta': DIHEDRAL_BINS}
    #
    hist = {}
    for feature_name in feature_s:
        hist[feature_name] = np.zeros((n_residue, n_residue, BINS[feature_name].shape[0]), dtype=np.float32)
        #
        for i in range(n_residue):
            for j in range(n_residue):
                if i == j: continue
                #
                if feature_name in ['dist', 'omega'] and j < i:
                    hist[feature_name][j,i] = hist[feature_name][i,j]
                    continue
                #
                h = np.histogram(feature_s[feature_name][:,i,j], bins=BINS[feature_name])[0].astype(np.float32)
                h /= h.sum()
                hist[feature_name][i,j,0] = h[-1]
                hist[feature_name][i,j,1:] = h[:-1]
    return hist

def gaussian(m, s, X0):
    return np.exp(-0.5*((X0[None,:]-m[:,None])/s)**2) / (np.sqrt(2.*np.pi) * s)

def vonMises(m, kappa, X0):
    return np.exp(kappa*np.cos(X0[None,:]-m[:,None])) / (2.*np.pi*np.i0(kappa))

def get_distr(n_residue, feature_s, mask):
    TYPE = {'dist': 'distance', 'phi': 'angle', 'omega': 'dihedral', 'theta': 'dihedral'}
    SIGMA = {'distance': 0.5, 'angle': 20., 'dihedral': 20.}
    SPACE = {'distance': 0.5, 'angle': np.deg2rad(15.), 'dihedral': np.deg2rad(15.)}
    #
    DIST_BINS = np.linspace(2., 20., 37) + SPACE['distance']*0.5
    ANGLE_BINS = np.deg2rad(np.linspace(0., 180., 13)) + SPACE['angle']*0.5
    DIHEDRAL_BINS = np.deg2rad(np.linspace(-180., 180., 25)) + SPACE['dihedral']*0.5
    BINS = {'distance': DIST_BINS, 'angle': ANGLE_BINS, 'dihedral': DIHEDRAL_BINS}
    FUNCTIONS = {'distance': gaussian, 'angle': vonMises, 'dihedral': vonMises}
    #
    hist = {}
    for feature_name in feature_s:
        feature_type = TYPE[feature_name]
        space = SPACE[feature_type]
        sigma = SIGMA[feature_type]
        bins = BINS[feature_type]
        func = FUNCTIONS[feature_type]
        #
        hist[feature_name] = np.zeros((n_residue, n_residue, bins.shape[0]), dtype=np.float32)
        #
        for i in range(n_residue):
            if not mask[i]: continue
            for j in range(n_residue):
                if not mask[j]: continue
                if i == j: continue
                #
                if feature_name in ['dist', 'omega'] and j < i:
                    hist[feature_name][i,j] = hist[feature_name][j,i]
                    continue
                #
                m = feature_s[feature_name][:,i,j]
                h0 = func(m, sigma, bins) * space
                h0 = h0.sum(axis=0)
                h0 /= np.float32(m.shape[0])
                h = np.zeros_like(h0)
                h[1:] = h0[:-1]
                h[0] = 1.-h[1:].sum()
                h = np.maximum(0, h)
                hist[feature_name][i,j] = h
    return hist

def read_summary(fn):
    prefix_s = []
    with fn.open() as fp:
        for line in fp:
            x = line.strip().split()
            prefix_s.append(x[4])
    return prefix_s

def read_pir(fn):
    query = []
    templ_s = []
    with fn.open() as fp:
        is_query = None
        for line in fp:
            if line.startswith(">"): 
                is_query = None
                continue
            if line.startswith("sequence:"):
                is_query = True
            elif line.startswith("structure:"):
                is_query = False
                templ = [] ; templ_s.append(templ)
            elif is_query is not None:
                seq = line.strip().replace("*","")
                if is_query:
                    query.append(seq)
                else:
                    templ.append(seq)
    #
    query = ''.join(query)
    templ_s = [''.join(templ) for templ in templ_s]
    #
    query_mask = (np.array(list(query), dtype='<U1') != '-')
    templ_mask_s = np.array([(np.array(list(templ), dtype='<U1') != '-') for templ in templ_s])
    #
    mask = np.any(templ_mask_s, axis=0)[query_mask]

    return mask

def hybrid_contacts(weight_s, distr_s):
    hybrid = {}
    for feature in FEATUREs:
        hybrid[feature] = np.zeros_like(distr_s[0][feature])
        l_seq = hybrid[feature][0].shape[0]
        weight_sum = np.zeros((l_seq, l_seq), dtype=np.float32)
        #
        for weight, distr in zip(weight_s, distr_s):
            w = np.ix_(weight, weight)
            weight_sum[w] += 1.0
            hybrid[feature] += distr[feature]
        hybrid[feature] /= weight_sum[:,:,None]
        zero = (weight_sum == 0.0)
        hybrid[feature][zero] = 0.0
    return hybrid, weight_sum

def run(title, domain):
    WEIGHTs = {'blast': 0.4, 'hhpred': 0.3, 'hmmer': 0.3}
    #
    out_s = []
    for method in METHODs:
        method_home = domain.domain_home().subdir(method)
        summary = read_summary(method_home.fn("%s.templ_s.summary"%(title)))
        weight = 0.1*WEIGHTs[method]
        if len(summary) == 0:
            continue
        #
        out_fn = method_home.fn("%s.tbm3.npz"%(title))
        if out_fn.status():
            out = np.load(out_fn.short())
            weight_sum = out['weight']
            out_s.append((weight_sum, out))
            continue
        #
        weight_per_residue = []
        distr_s = []
        for prefix in summary:
            pdb_fn_s = method_home.glob("%s.B*.pdb"%(prefix))
            if len(pdb_fn_s) == 0: continue
            #
            pir_fn = method_home.fn("%s.pir"%prefix)
            mask = read_pir(pir_fn)
            weight_per_residue.append(mask)
            #
            pdb_s = []
            for pdb_fn in pdb_fn_s:
                pdb = mdtraj.load(pdb_fn.short())
                pdb = pdb.atom_slice(pdb.top.select("name CA or name CB or name N or name O"))
                pdb_s.append(pdb)
            pdb_s = mdtraj.join(pdb_s)
            #
            n_residue = pdb_s.top.n_residues
            is_gly, feature_s = get_features(pdb_s)
            distr_s.append(get_distr(n_residue, feature_s, mask))
        if len(distr_s) == 0: continue
        #
        out, weight_sum = hybrid_contacts(weight_per_residue, distr_s)
        weight_sum *= weight
        out['weight'] = weight_sum
        #
        np.savez(out_fn.short(), **out)
        out_s.append((weight_sum, out))
    return out_s

