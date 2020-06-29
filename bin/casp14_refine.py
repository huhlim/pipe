#!/usr/bin/env python

import os
import sys
import path
import argparse
import subprocess as sp
from importlib import import_module

from libcommon import *
from libmain import *
from libligand import get_ligand_info

def get_membrane_topology(job, n_init, wait_after_run, sleep=30):
    membrane_home = job.work_home.subdir("membrane", build=True)
    #
    while True:
        job.membrane_pdb = []
        job.membrane_psf = []
        job.membrane_crd = []
        #
        status = True
        for i in range(n_init):
            m_home = membrane_home.subdir("%d"%i, build=True)
            pdb_fn = m_home.glob("*.pdb")
            psf_fn = m_home.glob("*.psf")
            crd_fn = m_home.glob("*.crd")
            if len(pdb_fn) == 0:
                status = False ; break
            if len(psf_fn) == 0:
                status = False ; break
            if len(crd_fn) == 0:
                status = False ; break
            #
            job.membrane_pdb.append(pdb_fn[0])
            job.membrane_psf.append(psf_fn[0])
            job.membrane_crd.append(crd_fn[0])
        #
        if wait_after_run and (not status):
            sys.stderr.write("waiting for CHARMM-GUI membrane topology... \n")
            time.sleep(sleep)
        else:
            break
        #
    if status: job.to_json()
    return status

def get_oligomer_topology(job, wait_after_run, sleep=30):
    oligomer_home = job.work_home.subdir("oligomer", build=True)
    #
    while True:
        oligomer_pdb = oligomer_home.glob("*.pdb")
        #
        status = (len(oligomer_pdb) != 0)
        #
        if wait_after_run and (not status):
            sys.stderr.write("waiting for Oligomer topology... \n")
            time.sleep(sleep)
        else:
            break
        #
    if status: 
        job.oligomer_pdb = oligomer_pdb
        job.to_json()
    return status

def main():
    arg = argparse.ArgumentParser(prog='PREFMD')
    arg.add_argument(dest='title', help='Job title')
    arg.add_argument('-i', '--input', dest='input_pdb', \
            help='input PDB file')
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
    arg.add_argument('--extensive', dest='use_extensive', action='store_true', default=False, \
            help='use extensive sampling')
    arg.add_argument('--membrane', dest='is_membrane_protein', action='store_true', default=False)
    arg.add_argument('--ligand', dest='has_ligand', action='store_true', default=False)
    arg.add_argument('--oligomer', dest='is_oligomer', action='store_true', default=False)

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    if arg.input_pdb is not None:
        arg.input_pdb = path.Path(arg.input_pdb)
    #
    # init
    if arg.title.endswith("job.json"):
        input_json = path.Path(arg.title)
        job = Job.from_json(input_json)
        job.verbose = arg.verbose
        job.keep_tmp = arg.keep
        job.append_to_joblist()
    else:
        job = import_module("init_refine").prep(arg)
    #
    if job.has("is_oligomer") and not job.has("oligomer_pdb"):
        if not get_oligomer_topology(job, arg.wait_after_run):
            sys.stderr.write("waiting for Oligomer topology... \n")
            return

    # hybrid
    if job.use_hybrid:
        import_module("hybrid").prep(job, job.init_pdb[0])
        if not run(job, arg.wait_after_run):
            return 
        hybrid_out = get_outputs(job, 'hybrid', expand='model_s')
        job.init_pdb.extend(hybrid_out[0][0])
        job.init_pdb = job.init_pdb[:5]
        job.to_json()
    #
    n_init = len(job.init_pdb)

    # locPREFMD
    if job.has("is_oligomer"):
        import_module("locPREFMD").prep(job, job.oligomer_pdb)
    else:
        import_module("locPREFMD").prep(job, job.init_pdb)
    if not run(job, arg.wait_after_run):
        return 
    locPREFMD_out = get_outputs(job, "locPREFMD")[:n_init]

    # define topology
    import_module("define_topology").prep(job, locPREFMD_out[0][0])
    if job.has("is_membrane_protein"):
        if not get_membrane_topology(job, n_init, arg.wait_after_run):
            sys.stderr.write("waiting for CHARMM-GUI membrane topology... \n")
            return
    if job.has("has_ligand"):
        if not get_ligand_info(job, arg.wait_after_run):
            sys.stderr.write("waiting for ligand info... \n")
            return

    # equil
    if not job.has("is_membrane_protein"):
        import_module("equil").prep(job, 0, [out[0] for out in locPREFMD_out], \
                                                path.Path("%s/equil.json"%(DEFAULT_HOME)))
    else:
        import_module("equil").prep_membrane(job, 0, \
                                            job.membrane_pdb, job.membrane_psf, job.membrane_crd, \
                                            path.Path("%s/equil_membrane.json"%(DEFAULT_HOME)))
    if not run(job, arg.wait_after_run):
        return 
    
    # prod
    if not job.has("is_membrane_protein"):
        if job.use_extensive:
            prod_input = path.Path("%s/prod_ext.json"%DEFAULT_HOME)
            n_traj = 20
        else:
            prod_input = path.Path("%s/prod.json"%DEFAULT_HOME)
            n_traj = 5
        for i in range(n_init):
            import_module("prod").prep(job, i, i, prod_input, n_traj)
    else:
        prod_input = path.Path("%s/prod_membrane.json"%DEFAULT_HOME)
        n_traj = 5
        for i in range(n_init):
            import_module("prod").prep(job, i, i, prod_input, n_traj)
    if not run(job, arg.wait_after_run):
        return
    prod_out = get_outputs(job, 'prod')
    
    # score
    import_module("score").prep(job, [out[0] for out in prod_out])
    if not run(job, arg.wait_after_run):
        return

    # average
    if n_init == 1:
        import_module("average").prep(job, '%s.prod_0'%job.title, [0], path.Path("%s/average.json"%DEFAULT_HOME), rule='score')
        import_module("average").prep(job, '%s.prod_0.cluster'%job.title, [0], path.Path("%s/average.json"%DEFAULT_HOME), rule='cluster')
        if not run(job, arg.wait_after_run):
            return
    else:
        import_module("average").prep(job, job.title, [i for i in range(n_init)], path.Path("%s/average.json"%DEFAULT_HOME), rule='casp12')
        for i in range(n_init):
            import_module("average").prep(job, '%s.prod_%d'%(job.title, i), [i], path.Path("%s/average.json"%DEFAULT_HOME), rule='score')
        if not run(job, arg.wait_after_run):
            return
    average_out = []
    for output in get_outputs(job, 'average', expand='pdb_s'):
        average_out.extend(output[0])
    average_out = average_out[:N_MODEL]

    # model
    job.work_home.chdir()
    model_home = job.work_home.subdir("model", build=True)
    prep_s = []
    for i,out in enumerate(average_out):
        prep_fn = model_home.fn("prep_%d.pdb"%(i+1))
        if not prep_fn.status():
            system(['cp', out.short(), prep_fn.short()], verbose=job.verbose)
        prep_s.append(prep_fn)
    #
    import_module("scwrl").prep(job, prep_s)
    if not run(job, arg.wait_after_run):
        return
    scwrl_out = get_outputs(job, 'scwrl')
    #
    import_module("locPREFMD").prep(job, [out[0] for out in scwrl_out])
    if not run(job, arg.wait_after_run):
        return 
    locPREFMD_out = get_outputs(job, "locPREFMD")[n_init:]
    #
    model_s = []
    for i,out in enumerate(locPREFMD_out):
        model_fn = model_home.fn("model_%d.pdb"%(i+1))
        if not model_fn.status():
            system(['cp', out[0].short(), model_fn.short()], verbose=job.verbose)
        model_s.append(model_fn)

    # qa
    import_module("qa").prep(job, model_s, path.Path("%s/qa.json"%DEFAULT_HOME))
    if not run(job, arg.wait_after_run):
        return
    qa_out = get_outputs(job, 'qa')

    # final
    final_home = job.work_home.subdir("final", build=True)
    for i,out in enumerate(qa_out):
        pdb_fn = final_home.fn("model_%d.pdb"%(i+1))
        if not pdb_fn.status():
            system(['cp', out[0].short(), pdb_fn.short()], verbose=job.verbose)
    #
    job.remove_from_joblist()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit()

