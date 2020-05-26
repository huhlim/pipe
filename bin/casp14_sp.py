#!/usr/bin/env python

import os
import sys
import path
import argparse
import numpy as np
import subprocess as sp
from importlib import import_module

from libcommon import *
from libmain import *

EXEC_REFINE = '%s/casp14_refine.py'%BIN_HOME
EXEC_PASTE = '%s/trRosetta/paste_domains.py'%EXEC_HOME

N_MODEL_REFINE = 5
N_MODEL = 5

PARAM_HYBRID_REFINE_MAX_RES = 300

def run_refine(title, input_pdb, work_home, **kwargs):
    def get_num_residues(pdb_fn):
        n_res = 0
        with pdb_fn.open() as fp:
            for line in fp:
                if not line.startswith("ATOM"):
                    continue
                atmName = line[12:16].strip()
                if atmName == 'CA':
                    n_res += 1
        return n_res
    def is_continuous_domain(pdb_fn):
        R = []
        with pdb_fn.open() as fp:
            for line in fp:
                if not line.startswith("ATOM"):
                    continue
                atmName = line[12:16].strip()
                if atmName == 'CA':
                    R.append([line[30:38], line[38:46], line[46:54]])
        R = np.array(R, dtype=float)
        b = np.sqrt(np.sum((R[1:] - R[:-1])**2, -1))
        b_max = np.max(b)
        return b_max < 5.0
    #
    n_residue = get_num_residues(input_pdb)
    #
    cmd = []
    cmd.append(EXEC_REFINE)
    cmd.append(title)
    cmd.extend(['--input', input_pdb.short()])
    cmd.extend(['--dir', work_home.short()])
    if kwargs.get("verbose", False):
        cmd.append("--verbose")
    if kwargs.get("wait_after_run", False):
        cmd.append("--wait")
    if kwargs.get("is_membrane_protein", False):
        cmd.append("--membrane")
    if kwargs.get("has_ligand", False):
        cmd.append("--ligand")
    if kwargs.get("use_hybrid", False) and \
            is_continuous_domain(input_pdb) and \
            (n_residue < PARAM_HYBRID_REFINE_MAX_RES):
        cmd.append("--hybrid")
    #
    proc = sp.Popen(cmd)
    time.sleep(60)
    return proc

def wait_refine(refine_proc_s, refine_home_s, wait_after_run, sleep=30):
    # check proc_s
    if refine_proc_s is not None:
        while True:
            proc_status = True
            for refine_proc in refine_proc_s:
                if refine_proc.poll() is None:
                    proc_status = False
            #
            if proc_status:
                break
            #
            time.sleep(sleep)
    #
    # check model_s
    status = True
    refined = []
    for refine_home in refine_home_s:
        model_home = refine_home.subdir("final")
        model_s = model_home.glob("model_?.pdb")
        model_s.sort(key=lambda x: int(x.name().split("_")[1]))
        refined.append(model_s)
        if len(model_s) < N_MODEL_REFINE:
            status = False
    #
    if wait_after_run and status:
        # expect to have refined models after all the process ends
        # , so if not raise error.
        sys.exit("Error: failed to get refined models\n")

    return status, refined

def paste_refined(model_home, framework, refined_s, verbose=False):
    product = import_module("itertools").product
    #
    n_model = [np.arange(len(X), dtype=int) for X in refined_s]
    product_s = [X for X in product(*n_model)]
    product_s.sort(key=lambda x: sum(x))
    #
    prep_s = []
    for i in range(N_MODEL):
        prep_fn = model_home.fn("prep_%d.pdb"%(i+1))
        prep_s.append(prep_fn)
        if prep_fn.status():
            continue
        #
        cmd = [EXEC_PASTE]
        cmd.append("--init")
        cmd.append(framework.short())
        cmd.append("--refined")
        for k,refined in zip(product_s[i], refined_s):
            cmd.append(refined[k].short())
        cmd.append("--output")
        cmd.append(prep_fn.short())
        system(cmd, verbose=verbose)
        if not prep_fn.status():
            sys.exit("ERROR: Failed to paste models\n")
    return prep_s

def update_bfactor(out_fn, model_fn, bfac_fn):
    bfac_s = {}
    with bfac_fn.open() as fp:
        for line in fp:
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                continue
            atmName = line[12:16].strip()
            if atmName != 'CA':
                continue
            resNo = int(line[22:26])
            bfac = float(line[60:66])
            bfac_s[resNo] = bfac
    #
    wrt = []
    with model_fn.open() as fp:
        for line in fp:
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                wrt.append(line)
                continue
            resNo = int(line[22:26])
            bfac = bfac_s[resNo]
            wrt.append('%s%6.2f\n'%(line[:60], bfac))
    #
    with out_fn.open('wt') as fout:
        fout.writelines(wrt)

def main():
    arg = argparse.ArgumentParser(prog='trRosetta+PREFMD')
    arg.add_argument(dest='title', help='Job title')
    arg.add_argument('-i', '--input', dest='input_fa', \
            help='input FA file')
    arg.add_argument('-d', '--dir', dest='work_dir', default='./',\
            help='working directory (default=./)')
    arg.add_argument('--keep', dest='keep', action='store_true', default=False,\
            help='set temporary file mode (default=False)')
    arg.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,\
            help='set verbose mode (default=False)')
    arg.add_argument('-w', '--wait', dest='wait_after_run', action='store_true', default=False,\
            help='set running type (default=False)')
    arg.add_argument('--hybrid', dest='use_hybrid', action='store_true', default=False, \
            help='use hybrid')
    arg.add_argument('--membrane', dest='is_membrane_protein', action='store_true', default=False)
    arg.add_argument('--ligand', dest='has_ligand', action='store_true', default=False)

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    if arg.input_fa is not None:
        arg.input_fa = path.Path(arg.input_fa)
    #
    # init
    if arg.title.endswith("job.json"):
        input_json = path.Path(arg.title)
        job = Job.from_json(input_json)
        job.verbose = arg.verbose
        job.keep_tmp = arg.keep
        job.append_to_joblist()
    else:
        job = import_module("init_sp").prep(arg)
    
    # trRosetta
    import_module("trRosetta").prep(job, job.init_fa)
    if not run(job, arg.wait_after_run):
        return 
    domain_pdb_s, trRosetta_min = get_outputs(job, 'trRosetta', expand='model_s')[0]

    # create refine directory
    job.refine_home = job.work_home.subdir("refine", build=True)
    #
    has_refine = job.has("refine_s")
    if has_refine:
        for refine_home in job.refine_s:
            if not refine_home.fn("job.json").status():
                has_refine = False ; break
    #
    if not has_refine:
        refine_proc_s = []
        job.refine_s = []
        for pdb_fn in domain_pdb_s:
            domain_id = pdb_fn.name()
            #
            refine_proc = run_refine(domain_id, pdb_fn, job.refine_home, verbose=arg.verbose, \
                                     use_hybrid=job.use_hybrid, \
                                     is_membrane_protein=job.has("is_membrane_protein"), \
                                     has_ligand=job.has("has_ligand"), \
                                     wait_after_run=arg.wait_after_run)
            refine_home = job.refine_home.subdir(domain_id)
            #
            refine_proc_s.append(refine_proc)
            job.refine_s.append(refine_home)
    else:
        refine_proc_s = None
    job.to_json()
    #
    status, refined = wait_refine(refine_proc_s, job.refine_s, arg.wait_after_run)
    if not status:
        return
    #
    # paste
    job.work_home.chdir()
    model_home = job.work_home.subdir("model", build=True)
    prep_s = paste_refined(model_home, trRosetta_min, refined, verbose=arg.verbose)
    #
    import_module("scwrl").prep(job, prep_s)
    if not run(job, arg.wait_after_run):
        return
    scwrl_out = get_outputs(job, 'scwrl')
    #
    import_module("locPREFMD").prep(job, [out[0] for out in scwrl_out])
    if not run(job, arg.wait_after_run):
        return 
    locPREFMD_out = get_outputs(job, "locPREFMD")

    # final
    final_home = job.work_home.subdir("final", build=True)
    for i,out in enumerate(locPREFMD_out):
        pdb_fn = final_home.fn("model_%d.pdb"%(i+1))
        if not pdb_fn.status():
            update_bfactor(pdb_fn, out[0], prep_s[i])
    #
    job.remove_from_joblist()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit()

