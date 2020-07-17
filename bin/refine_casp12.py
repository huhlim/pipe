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
    arg.add_argument('--ligand', dest='has_ligand', action='store_true', default=False)

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    if arg.input_pdb is not None:
        arg.input_pdb = path.Path(arg.input_pdb)
    #
    arg.use_hybrid = False
    arg.use_extensive = False
    arg.is_membrane_protein = False
    arg.is_oligomer = False
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
        job.run_type = 'refine_casp12'
        job.run_exec=os.path.abspath(__file__)
    #
    n_init = len(job.init_pdb)

    # locPREFMD
    import_module("locPREFMD").prep(job, job.init_pdb)
    if not run(job, arg.wait_after_run):
        return 
    locPREFMD_out = get_outputs(job, "locPREFMD")[:n_init]

    # define topology
    import_module("define_topology").prep(job, locPREFMD_out[0][0])
    if job.has("has_ligand"):
        if not get_ligand_info(job, arg.wait_after_run):
            sys.stderr.write("waiting for ligand info... \n")
            return

    # equil
    import_module("equil").prep(job, 0, [out[0] for out in locPREFMD_out], \
                                            path.Path("%s/equil_casp12.json"%(DEFAULT_HOME)))
    if not run(job, arg.wait_after_run):
        return 
    
    # prod
    n_traj = 5
    prod_input = path.Path("%s/prod_casp12.json"%DEFAULT_HOME)
    for i in range(n_init):
        import_module("prod").prep(job, i, i, prod_input, n_traj)
    if not run(job, arg.wait_after_run):
        return
    prod_out = get_outputs(job, 'prod')
    #
    job.remove_from_joblist()
    return
    #
    # score
    import_module("score").prep(job, [out[0] for out in prod_out])
    if not run(job, arg.wait_after_run):
        return

    # average
    import_module("average").prep(job, '%s.casp12'%job.title,  [0], path.Path("%s/average.json"%DEFAULT_HOME), rule='casp12')
    import_module("average").prep(job, '%s.score'%job.title,   [0], path.Path("%s/average.json"%DEFAULT_HOME), rule='score')
    import_module("average").prep(job, '%s.cluster'%job.title, [0], path.Path("%s/average.json"%DEFAULT_HOME), rule='cluster')
    if not run(job, arg.wait_after_run):
        return
    #
    average_out = []
    for output in get_outputs(job, 'average', expand='pdb_s'):
        average_out.extend(output[0])
    #average_out = average_out[:N_MODEL]

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

    # final
    final_home = job.work_home.subdir("final", build=True)
    for i,out in enumerate(locPREFMD_out):
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

