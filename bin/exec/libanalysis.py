#!/usr/bin/env python

import os
import sys
import path
import mdtraj
import numpy as np
import pickle

def get_molecules(topo, exclude_s=[], max_water=-1):
    if isinstance(topo, mdtraj.Trajectory):
        top = topo.topology
    elif isinstance(topo, mdtraj.Topology):
        top = topo
    else:
        raise ValueError(topo)
    #
    molecule_s = {}
    for molecule in top.find_molecules():
        index = np.array(sorted([atom.index for atom in molecule]), dtype=int)
        residue = list(molecule)[0].residue
        if residue.is_protein:    # they have unique segment_ids
            name = residue.segment_id
        elif residue.is_water:
            name = 'water'
        else:   # ions
            name = residue.name
        if name in exclude_s:
            continue
        #
        if name not in molecule_s:
            molecule_s[name] = []
        molecule_s[name].append(index)
    if max_water >= 0 and 'water' in molecule_s:
        water = np.array(molecule_s['water'])
        index = np.random.choice(list(range(len(water))), size=max_water)
        molecule_s['water'] = water[index]

    return molecule_s

def calc_center_of_mass(molecule_s, traj, include_s=[], exclude_s=[]):
    if not isinstance(traj, mdtraj.Trajectory):
        raise ValueError(traj)
    #
    if len(include_s) == 0:
        molecule_name_s = list(molecule_s.keys())
    else:
        molecule_name_s = []
        for name in include_s:
            if name in molecule_s:
                molecule_name_s.append(name)
            else:
                raise KeyError(name)
    molecule_name_s = [name for name in molecule_name_s if name not in exclude_s]
    #
    mass = np.array([atom.element.mass for atom in traj.top.atoms], dtype=float)
    #
    center_of_mass = {}
    for molecule_name in molecule_name_s:
        index_s = molecule_s[molecule_name]
        cntr_s = []
        for index in index_s:
            xyz = traj.xyz[:,index]
            m = mass[index]
            cntr = (xyz*m[None,:,None]).sum(axis=1) / m.sum()
            cntr_s.append(cntr)
        #
        center_of_mass[molecule_name] = np.array(cntr_s, dtype=float)   # shape=(# molecule, # frames, 3)
    #
    return center_of_mass

def calc_average_boxsize(traj_s):
    pbc_s = []
    for traj in traj_s:
        if isinstance(traj, mdtraj.Trajectory):
            pbc = traj.unitcell_lengths
        else:
            pbc = np.concatenate([t.unitcell_lengths for t in traj])
        pdb_s.append(np.mean(pbc, axis=0))
    #
    mean = np.mean(pbc_s, axis=0)
    std = np.std(pbc_s, axis=0)
    return mean, std

def unwrap_PBC(r, pbc0):
    pbc = pbc0[1:][None,:]
    dr = r[:,1:] - r[:,:-1]
    for i in range(3):
        wrap = np.cumsum(dr[:,:,i] <-pbc[:,:,i]*0.5, axis=1) - np.cumsum(dr[:,:,i] > pbc[:,:,i]*0.5, axis=1)
        r[:,1:,i] += wrap*pbc[:,:,i]
    return r
