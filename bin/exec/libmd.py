#!/usr/bin/env python

import os
import sys

import warnings
warnings.filterwarnings("ignore")

import numpy as np
from sklearn.decomposition import PCA

import mdtraj
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

WORK_HOME = os.getenv("PREFMD_HOME")
assert WORK_HOME is not None
sys.path.insert(0, '%s/bin'%WORK_HOME)

import path
from libcommon import *

from libquat import Quaternion

def solvate_pdb(output_prefix, pdb, options, verbose):
    orient_fn = path.Path('%s.orient.pdb'%output_prefix)
    if not orient_fn.status():
        xyz = pdb.xyz[0]
        xyz -= xyz.mean(axis=0)
        #
        for i in range(2):
            pca = PCA(n_components=(i+1))
            pca.fit(xyz)
            #
            axis_0 = pca.components_[i]
            axis_1 = np.zeros(3, dtype=float)
            axis_1[i] = 1.
            #
            axis_r = np.cross(axis_0, axis_1)
            angle_r = np.arccos(np.dot(axis_0, axis_1))
            #
            q = Quaternion.from_axis_and_angle(axis_r, angle_r)
            #
            xyz = np.dot(xyz, q.rotate().T)
        #
        system_size = np.max(xyz, axis=0) - np.min(xyz, axis=0) + 2.0*(options['md']['solvate'] * 0.1)
        translate = (system_size - (np.max(xyz, axis=0) + np.min(xyz, axis=0))) * 0.5
        xyz += translate
        pdb.xyz[0] = xyz
        #
        pdb.unitcell_vectors = (system_size * np.eye(3))[None,:]
        pdb.save(orient_fn.short())
    #
    solv_fn = path.Path('%s.solvate.pdb'%output_prefix)
    if not solv_fn.status():
        cmd = []
        cmd.append("%s/solvate.py"%EXEC_HOME)
        cmd.append(orient_fn.short())
        cmd.append(solv_fn.short())
        cmd.append("%8.5f"%options['md']['ion_conc'])
        system(cmd, verbose=verbose)
        #
        cmd = []
        cmd.append("%s/update_water_name.py"%EXEC_HOME)
        cmd.append(solv_fn.short())
        system(cmd, verbose=verbose)
        #
        if 'use_modified_CMAP' in options['ff'] and options['ff']['use_modified_CMAP']:
            cmd = ['%s/resName_modified_CMAP.py'%EXEC_HOME, solv_fn.short()]
            output = system(cmd, verbose=verbose, stdout=True)
        else:
            with solv_fn.open() as fp:
                output = fp.read()
        #
        with solv_fn.open('wt') as fout:
            if 'ssbond' in options:
                for line in options['ssbond']:
                    fout.write("%s\n"%line)
            fout.write(output)
    #
    return orient_fn, solv_fn

def generate_PSF(output_prefix, solv_fn, options, verbose):
    psf_fn = path.Path("%s.psf"%output_prefix)
    crd_fn = path.Path("%s.crd"%output_prefix)
    if psf_fn.status() and crd_fn.status():
        return psf_fn, crd_fn
    #
    cmd = []
    cmd.append("%s/genPSF.py"%EXEC_HOME)
    cmd.append(solv_fn.short())
    cmd.extend(['-psf', psf_fn.short()])
    cmd.extend(['-crd', crd_fn.short()])
    cmd.append("--toppar")
    cmd.extend(options['ff']['toppar'])
    if 'ligand' in options and len(options['ligand']['str_fn_s']) > 0:
        cmd.extend(options['ff']['cgenff'])
        cmd.extend([fn.short() for fn in options['ligand']['str_fn_s']])
    system(cmd, verbose=verbose)
    #
    return psf_fn, crd_fn

def construct_restraint(psf, pdb, force_const):
    rsr = CustomExternalForce("k0*d^2 ; d=periodicdistance(x,y,z, x0,y0,z0)")
    rsr.addPerParticleParameter("x0")
    rsr.addPerParticleParameter("y0")
    rsr.addPerParticleParameter("z0")
    rsr.addPerParticleParameter('k0')
    #
    calphaIndex = []
    for i,atom in enumerate(psf.topology.atoms()):
        if atom.name == 'CA':
            calphaIndex.append(i)
    #
    k = -1
    for i,atom in enumerate(pdb.top.atoms):
        if atom.name != 'CA': continue
        #
        k += 1
        mass = atom.element.mass
        param = pdb.xyz[0,i].tolist()
        param.append(force_const * mass * kilocalories_per_mole/angstroms**2)
        rsr.addParticle(calphaIndex[k], param)
    return rsr

def construct_membrane_restraint(psf, pdb, force_const):
    rsr = CustomExternalForce("k0*d^2 ; d=periodicdistance(x,y,z, x0,y0,z0)")
    rsr.addPerParticleParameter("x0")
    rsr.addPerParticleParameter("y0")
    rsr.addPerParticleParameter("z0")
    rsr.addPerParticleParameter('k0')
    #
    heavyIndex = []
    for i,atom in enumerate(psf.topology.atoms()):
        #if atom.element.mass > 4.0*amu:
        if atom.name == 'P':
            heavyIndex.append(i)
    #
    k = -1
    for i,atom in enumerate(pdb.top.atoms):
        #if atom.element.mass < 4.0:
        if atom.name != 'P':
            continue
        #
        k += 1
        mass = atom.element.mass
        param = pdb.xyz[0,i].tolist()
        param.append(force_const * mass * kilocalories_per_mole/angstroms**2)
        rsr.addParticle(heavyIndex[k], param)
    return rsr

def construct_ligand_restraint(pair_s):
    bond = CustomBondForce("k * (r-r0)^2")
    bond.addPerBondParameter('k')
    bond.addPerBondParameter('r0')
    #
    for pair in pair_s:
        bond.addBond(pair[0], pair[1], \
                (pair[2]*kilocalories_per_mole/angstroms**2, pair[3]*nanometers))
    return bond
