#!/usr/bin/env python

import os
import sys
import path
import argparse
import subprocess as sp
from importlib import import_module

from libcommon import *
from libmain import *

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

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    if arg.input_pdb is not None:
        arg.input_pdb = path.Path(arg.input_pdb)
    #
    # init
    if arg.title.endswith("/job.json"):
        input_json = path.Path(arg.title)
        job = Job.from_json(input_json)
        job.verbose = arg.verbose
        job.keep_tmp = arg.keep
        job.append_to_joblist()
    else:
        job = import_module("init_refine").prep(arg)

    # hybrid
    if arg.use_hybrid:
        import_module("hybrid").prep(job, job.init_pdb[0])
        if not run(job, arg.wait_after_run):
            return 
        hybrid_out = get_outputs(job, 'hybrid')
        with hybrid_out[0][0].open() as fp:
            for line in fp:
                if line.startswith("#"): continue
                fn = path.Path(line.strip())
                job.init_pdb.append(fn)
        job.init_pdb = job.init_pdb[:5]
        job.to_json()
    #
    n_init = len(job.init_pdb)

    # locPREFMD
    import_module("locPREFMD").prep(job, job.init_pdb)
    if not run(job, arg.wait_after_run):
        return 
    locPREFMD_out = get_outputs(job, "locPREFMD")[:n_init]

    # define topology
    import_module("define_topology").prep(job, locPREFMD_out[0][0])

    # equil
    import_module("equil").prep(job, 0, [out[0] for out in locPREFMD_out], path.Path("%s/equil.json"%(DEFAULT_HOME)))
    if not run(job, arg.wait_after_run):
        return 
    
    # prod
    if arg.use_extensive:
        prod_input = path.Path("%s/prod_ext.json"%DEFAULT_HOME)
        n_traj = 10
    else:
        prod_input = path.Path("%s/prod.json"%DEFAULT_HOME)
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
    average_out = get_outputs(job, 'average', expand='pdb_s')[:N_MODEL]

    # model
    job.work_home.chdir()
    model_home = job.work_home.subdir("model", build=True)
    prep_s = []
    for i,out in enumerate(average_out):
        prep_fn = model_home.fn("prep_%d.pdb"%(i+1))
        if not prep_fn.status():
            system(['cp', out[0].short(), prep_fn.short()], verbose=arg.verbose)
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
            system(['cp', out[0].short(), model_fn.short()], verbose=arg.verbose)
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
            system(['cp', out[0].short(), pdb_fn.short()], verbose=arg.verbose)
    #
    job.remove_from_joblist()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit()

