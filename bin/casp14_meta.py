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

N_MODEL = 5

def add_meta_tasks(job, meta_s):
    method = 'prod'
    job.task[method] = []
    for meta in meta_s:
        job.task[method].extend([task for _,task in meta.get_task(method)])
    #
    method = 'score'
    job.task[method] = [] ; k = -1
    output_s = get_outputs(job, 'calc_rmsd')
    for i,meta in enumerate(meta_s):
        for j,task in meta.get_task(method):
            k += 1
            task['output'][1] = output_s[k][0]
            job.task[method].append(task)
    job.to_json()

def main():
    arg = argparse.ArgumentParser(prog='casp14_meta')
    arg.add_argument(dest='title', help='Job title')
    arg.add_argument('-i', '--input', dest='meta_json_s', nargs='*', \
            help='input JSON files', default=[])
    arg.add_argument('-d', '--dir', dest='work_dir', default='./',\
            help='working directory (default=./)')
    arg.add_argument('-w', '--wait', dest='wait_after_run', action='store_true', default=False,\
            help='set running type (default=False)')
    arg.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,\
            help='set verbose mode (default=False)')
    arg.add_argument('--keep', dest='keep', action='store_true', default=False,\
            help='set temporary file mode (default=False)')

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    arg.meta_json_s = [path.Path(fn) for fn in arg.meta_json_s]

    # init
    if arg.title.endswith("job.json"):
        input_json = path.Path(arg.title)
        job = Job.from_json(input_json)
        job.verbose = arg.verbose
        job.keep_tmp = arg.keep
        job.append_to_joblist()
    else:
        if len(arg.meta_json_s) == 0:
            sys.exit("Input json files were not given.\n")
        job = import_module("init_meta").prep(arg)
    #
    meta_s = []
    for meta_json_fn in job.meta_json_s:
        meta_s.append(Job.from_json(meta_json_fn))
    n_meta = len(meta_s)
    #
    job.init_pdb = [meta_s[0].init_pdb[0]]
    locPREFMD_out = get_outputs(meta_s[0], "locPREFMD")[0]
    import_module("define_topology").prep(job, locPREFMD_out[0])
    #
    # calc_rmsd
    for i,meta in enumerate(meta_s):
        prod_out = get_outputs(meta, 'prod')
        import_module("calc_rmsd").prep(job, i, job.init_pdb[0], [out[0] for out in prod_out])
    if not run(job, arg.wait_after_run, sleep=10):
        return

    add_meta_tasks(job, meta_s)

    # average
    import_module("average").prep(job, '%s.meta'%job.title, [i for i in range(n_meta)], path.Path("%s/average.json"%DEFAULT_HOME), rule='casp12')
    import_module("average").prep(job, '%s.cluster'%job.title, [i for i in range(n_meta)], path.Path("%s/average.json"%DEFAULT_HOME), rule='cluster')
    if not run(job, arg.wait_after_run):
        return
    #
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
            system(['cp', out.short(), prep_fn.short()], verbose=arg.verbose)
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

