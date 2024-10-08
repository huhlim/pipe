#!/usr/bin/env python

import os
import sys

import warnings

warnings.filterwarnings("ignore")

import numpy as np
from sklearn.decomposition import PCA

import mdtraj

from openmm.unit import *
from openmm.openmm import *
from openmm.app import *

try:
    from openmmtools.integrators import ThermostatedIntegrator
    from openmmtools.constants import kB

    openmmtools_available = True
except:
    openmmtools_available = False

WORK_HOME = os.getenv("PIPE_HOME")
assert WORK_HOME is not None
sys.path.insert(0, "%s/bin" % WORK_HOME)

import path
from libcommon import *

from libquat import Quaternion
from seqName import stdres_ext


def solvate_pdb(output_prefix, pdb, options, verbose):
    orient_fn = path.Path("%s.orient.pdb" % output_prefix)
    if not options["md"].get("orient", True):
        pdb.save(orient_fn.short())
    if not orient_fn.status():
        xyz = pdb.xyz[0]
        xyz -= xyz.mean(axis=0)
        #
        for i in range(2):
            pca = PCA(n_components=(i + 1))
            pca.fit(xyz)
            #
            axis_0 = pca.components_[i]
            axis_1 = np.zeros(3, dtype=float)
            axis_1[i] = 1.0
            if np.all(axis_0 == axis_1):
                continue
            #
            axis_r = np.cross(axis_0, axis_1)
            angle_r = np.arccos(np.dot(axis_0, axis_1))
            #
            q = Quaternion.from_axis_and_angle(axis_r, angle_r)
            #
            xyz = np.dot(xyz, q.rotate().T)
        #
        solvateBox = options["md"].get("solvateBox", "rectangular")
        if solvateBox == "cubic":
            system_size = np.ones(3) * max(np.max(xyz, axis=0) - np.min(xyz, axis=0))
        else:
            system_size = np.max(xyz, axis=0) - np.min(xyz, axis=0)
        system_size += 2.0 * (options["md"]["solvate"] * 0.1)
        translate = (system_size - (np.max(xyz, axis=0) + np.min(xyz, axis=0))) * 0.5
        xyz += translate
        pdb.xyz[0] = xyz
        #
        pdb.unitcell_vectors = (system_size * np.eye(3))[None, :]
        pdb.save(orient_fn.short())

    solvate_cryst = options["md"].get("solvate_cryst", None)
    if solvate_cryst is not None:
        cryst_pdb = path.Path("%s.cryst_water.pdb" % output_prefix)
        #
        if not cryst_pdb.status():
            cmd = []
            cmd.append("%s/copy_water.py" % EXEC_HOME)
            cmd.append(orient_fn.short())
            cmd.append(solvate_cryst[0].short())
            cmd.append(cryst_pdb.short())
            cmd.extend(["--bfactor", "%10.2f" % solvate_cryst[1]])
            system(cmd, verbose=verbose)
        #
        cryst_water = mdtraj.load(cryst_pdb.short())
        n_cryst_water = cryst_water.top.select("water").shape[0] // 3
        options["md"]["n_cryst_water"] = n_cryst_water
        #
        solv_input = cryst_pdb
    else:
        solv_input = orient_fn
    #
    solv_fn = path.Path("%s.solvate.pdb" % output_prefix)
    if not solv_fn.status():
        cmd = []
        cmd.append("%s/solvate.py" % EXEC_HOME)
        cmd.append(solv_input.short())
        cmd.append(solv_fn.short())
        cmd.append("%8.5f" % options["md"]["ion_conc"])
        system(cmd, verbose=verbose)
        #
        cmd = []
        cmd.append("%s/update_water_name.py" % EXEC_HOME)
        cmd.append(solv_fn.short())
        system(cmd, verbose=verbose)
        #
        if "use_modified_CMAP" in options["ff"] and options["ff"]["use_modified_CMAP"]:
            cmd = ["%s/resName_modified_CMAP.py" % EXEC_HOME, solv_fn.short()]
            output = system(cmd, verbose=verbose, stdout=True)
        else:
            with solv_fn.open() as fp:
                output = fp.read()
        #
        with solv_fn.open("wt") as fout:
            if "ssbond" in options:
                for line in options["ssbond"]:
                    fout.write("%s\n" % line)
            fout.write(output)
    #
    return orient_fn, solv_fn


def generate_PSF(output_prefix, solv_fn, options, verbose):
    psf_fn = path.Path("%s.psf" % output_prefix)
    crd_fn = path.Path("%s.crd" % output_prefix)
    if psf_fn.status() and crd_fn.status():
        return psf_fn, crd_fn
    #
    cmd = []
    cmd.append("%s/genPSF.py" % EXEC_HOME)
    cmd.append(solv_fn.short())
    cmd.extend(["-psf", psf_fn.short()])
    cmd.extend(["-crd", crd_fn.short()])
    cmd.append("--toppar")
    cmd.extend(options["ff"]["toppar"])
    if "blocked" in options["ff"] and options["ff"]["blocked"]:
        cmd.append("--blocked")
        if "terminal" in options["ff"]:
            cmd.append("--terminal")
            cmd.extend(options["ff"]["terminal"])
    if "patch" in options["ff"]:
        cmd.append("--patch")
        cmd.extend(options["ff"]["patch"])
    if "ligand" in options and len(options["ligand"]["str_fn_s"]) > 0:
        cmd.extend(options["ff"]["cgenff"])
        cmd.extend([fn.short() for fn in options["ligand"]["str_fn_s"]])
    system(cmd, verbose=verbose)
    #
    return psf_fn, crd_fn


def update_residue_name(pdb_fn, pdb):
    resName_s = ["HIS", "HSP", "HSD", "HSE"]
    residue_s = {}
    chain_s = []
    with pdb_fn.open() as fp:
        for line in fp:
            if not line.startswith("ATOM") and not line.startswith("HETA"):
                continue
            atmName = line[12:16].strip()
            if atmName != "CA":
                continue
            chain = line[21].strip()
            if chain not in chain_s:
                chain_s.append(chain)
            resName = line[17:20]
            if resName not in resName_s:
                continue
            resNo = line[22:27].strip()
            chain_index = chain_s.index(chain)
            try:
                residue_s[(chain_index, int(resNo))] = resName
            except:
                residue_s[(chain_index, resNo)] = resName
    #
    for residue in pdb.top.residues:
        chain_index = residue.chain.index
        resNo = residue.resSeq
        if (chain_index, resNo) in residue_s:
            resName = residue_s[(chain_index, resNo)]
            residue.name = resName


def put_segnames(in_pdb, CHAIN_BREAKs=2.0):
    seg_no = -1
    #
    segName_s = {}
    chain_prev = None
    R_prev = None
    with open(in_pdb) as fp:
        for line in fp:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            #
            resName = line[17:20]
            if resName not in stdres_ext:
                continue
            #
            atmName = line[12:16].strip()
            if atmName not in ["N", "C"]:
                continue
            #
            chain_id = line[21]
            resSeq = line[22:27].strip()
            key = (chain_id, resSeq)
            if key not in segName_s:
                segName_s[key] = None
                R = [None, None]
            #
            if chain_id != chain_prev:
                seg_no += 1
                xyz_prev = None
                chain_prev = chain_id
            #
            R[["N", "C"].index(atmName)] = np.array(
                [line[30:38], line[38:46], line[46:54]], dtype=float
            )
            if (R[0] is not None) and (R[1] is not None):
                xyz_curr = R[0]
                if xyz_prev is not None:
                    d = np.linalg.norm(xyz_curr - xyz_prev)
                    if d > CHAIN_BREAKs:
                        seg_no += 1
                xyz_prev = R[1]
                #
                segName_s[key] = f"P{seg_no:03d}"
    #
    out = []
    n_residue = {}
    key_prev = None
    with open(in_pdb) as fp:
        for line in fp:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                out.append(line)
                continue
            #
            line = line.rstrip()
            chain_id = line[21]
            resSeq = line[22:27].strip()
            key = (chain_id, resSeq)
            #
            resName = line[17:20]
            if resName in stdres_ext:
                segName = segName_s.get(key, None)
                if segName is None:
                    sys.exit(f"Failed to assign a segName to {chain_id}{resSeq}")
            else:
                if key != key_prev:
                    key_prev = key
                    if resName not in n_residue:
                        n_residue[resName] = 0
                    n_residue[resName] += 1
                    i_seg = n_residue[resName] // 10000

                if resName in ["SOD", "CLA", "POT"]:
                    segName = f"{resName}{i_seg:1d}"
                elif resName in ["TIP", "HOH"]:
                    segName = f"W{i_seg:03d}"
                else:
                    segName = f"{resName}{i_seg:1d}"
            #
            line = f"{line[:72]:<72s}{segName:<4s}{line[76:]}\n"
            out.append(line)

    return out, segName_s


def construct_restraint(psf, pdb, force_const, atom_s=["CA"]):
    rsr = CustomExternalForce("k0*d^2 ; d=periodicdistance(x,y,z, x0,y0,z0)")
    rsr.addPerParticleParameter("x0")
    rsr.addPerParticleParameter("y0")
    rsr.addPerParticleParameter("z0")
    rsr.addPerParticleParameter("k0")
    #
    calphaIndex = []
    for i, atom in enumerate(psf.topology.atoms()):
        if atom.name in atom_s:
            calphaIndex.append(i)
    #
    k = -1
    for i, atom in enumerate(pdb.top.atoms):
        if atom.name not in atom_s:
            continue
        #
        k += 1
        mass = atom.element.mass
        param = pdb.xyz[0, i].tolist()
        param.append(force_const * mass * kilocalories_per_mole / angstroms**2)
        rsr.addParticle(calphaIndex[k], param)
    return rsr


def construct_water_restraint(psf, pdb, n_water, force_const):
    rsr = CustomExternalForce("k0*d^2 ; d=periodicdistance(x,y,z, x0,y0,z0)")
    rsr.addPerParticleParameter("x0")
    rsr.addPerParticleParameter("y0")
    rsr.addPerParticleParameter("z0")
    rsr.addPerParticleParameter("k0")
    #
    waterIndex = []
    for residue in psf.topology.residues():
        if residue.name not in ["TIP3", "HOH"]:
            continue
        for atom in residue.atoms():
            if atom.name in ["OH2", "O"]:
                waterIndex.append(atom.index)
                break
        if len(waterIndex) >= n_water:
            break
    #
    pdb_water = pdb.atom_slice(pdb.top.select("water and element == O"))
    k = -1
    for i, atom in enumerate(pdb_water.top.atoms):
        k += 1
        mass = atom.element.mass
        param = pdb_water.xyz[0, i].tolist()
        param.append(force_const * mass * kilocalories_per_mole / angstroms**2)
        rsr.addParticle(waterIndex[k], param)
        if k + 1 == n_water:
            break
    return rsr


def construct_membrane_restraint(psf, pdb, force_const):
    rsr = CustomExternalForce("k0*d^2 ; d=periodicdistance(x,y,z, x0,y0,z0)")
    rsr.addPerParticleParameter("x0")
    rsr.addPerParticleParameter("y0")
    rsr.addPerParticleParameter("z0")
    rsr.addPerParticleParameter("k0")
    #
    heavyIndex = []
    for i, atom in enumerate(psf.topology.atoms()):
        # if atom.element.mass > 4.0*amu:
        if atom.name == "P":
            heavyIndex.append(i)
    #
    k = -1
    for i, atom in enumerate(pdb.top.atoms):
        # if atom.element.mass < 4.0:
        if atom.name != "P":
            continue
        #
        k += 1
        mass = atom.element.mass
        param = pdb.xyz[0, i].tolist()
        param.append(force_const * mass * kilocalories_per_mole / angstroms**2)
        rsr.addParticle(heavyIndex[k], param)
    return rsr


def construct_ligand_restraint(pair_s):
    bond = CustomBondForce("k * (r-r0)^2")
    bond.addPerBondParameter("k")
    bond.addPerBondParameter("r0")
    #
    for pair in pair_s:
        bond.addBond(
            pair[0],
            pair[1],
            (pair[2] * kilocalories_per_mole / angstroms**2, pair[3] * nanometers),
        )
    return bond


if openmmtools_available:

    class BerendsenVelocityVerletIntegrator(ThermostatedIntegrator):
        def __init__(
            self, temperature=298 * kelvin, timestep=1.0 * femtoseconds, risetime=1.0 * picoseconds
        ):
            super().__init__(temperature, timestep)

            self.addGlobalVariable("tau", timestep / risetime)
            self.addGlobalVariable("ke2", 0)  # Twice the kinetic energy
            self.addGlobalVariable("kTcurr", 0)
            self.addGlobalVariable("kB", kB)
            self.addGlobalVariable("lambda", 1.0)
            self.addGlobalVariable("ndf", 0)  # number of degrees of freedom
            self.addPerDofVariable("ones", 1.0)
            self.addPerDofVariable("x1", 0)
            ##
            self.addComputeSum("ndf", "ones")
            self.addComputeSum("ke2", "m*v*v")
            self.addComputeGlobal("kTcurr", "ke2/ndf*1.5")
            self.addComputeGlobal("lambda", "sqrt(1+tau*(kT/kTcurr - 1))")
            #
            self.addComputePerDof("v", "v*lambda")
            #
            # Velocity Verlet step
            #
            self.addUpdateContextState()
            self.addComputePerDof("v", "v+0.5*dt*f/m")
            self.addComputePerDof("x", "x+dt*v")
            self.addComputePerDof("x1", "x")
            self.addConstrainPositions()
            self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
            self.addConstrainVelocities()
            #


class PressureTensorReporter(object):
    def __init__(
        self, file, system, reportInterval, cacheInterval=None, enforcePeriodicBox=None, fmt="%.3f"
    ):
        self._reportInterval = reportInterval
        if cacheInterval is not None:
            self._cacheInterval = cacheInterval
        else:
            self._cacheInterval = reportInterval
        self._fmt = fmt
        self._enforcePeriodicBox = enforcePeriodicBox
        self._out = open(file, "wt")
        self._m = np.array(
            [
                system.getParticleMass(i).value_in_unit(dalton)
                for i in range(system.getNumParticles())
            ]
        )  # dalton = g/mol
        self._cache = []

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, True, True, False, self._enforcePeriodicBox)

    def report(self, simulation, state):
        F = state.getForces(asNumpy=True)  # kJ/mol/nm
        r = state.getPositions(asNumpy=True)  # nm
        v = state.getVelocities(asNumpy=True)  # nm/ps
        V = state.getPeriodicBoxVolume()._value
        #
        values = []
        values.append(str(simulation.currentStep))
        values.append(f"{state.getTime().value_in_unit(picosecond):.3f}")
        for (i, j) in [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]:
            p_ij_0 = np.dot(self._m * v[:, i], v[:, j])  # kJ/mol
            p_ij_1 = np.dot(r[:, i], F[:, j])  # kJ/mol
            p_ij = (p_ij_0 + p_ij_1) / V  # kJ/mol/nm^3
            values.append(self._fmt % p_ij)
        self._cache.append(values)
        if len(self._cache) >= self._cacheInterval:
            self.report_cache()

    def report_cache(self):
        if len(self._cache) >= self._cacheInterval:
            for values in self._cache:
                print(",".join(str(v) for v in values), file=self._out)
            try:
                self._out.flush()
            except AttributeError:
                pass
            self._cache = []

    def __del__(self):
        if len(self._cache) != 0:
            self.report_cache()
        self._out.close()
