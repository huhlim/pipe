#!/usr/bin/env python

import os
import sys
import path
import json
import time
import argparse
import subprocess as sp
from importlib import import_module

from libgpu import GPU
from libcommon import *

def check_status(job):
    updated = False
    status = True
    for method in job.task:
        task_s = job.get_task(method)
        for index, task in task_s:
            if task['resource'][0] == 'DONE':
                continue
            #
            task_done = True
            for output in task['output']:
                if not output.status():
                    task_done = False
                    break
            if task_done:
                updated = True
                job.update_task_status(method, index, 'DONE')
            else:
                status = False
    return updated, status

def get_resource_taken():
    cpu_s = [] ; gpu_s = {}
    if not JOBs_json.status():
        return cpu_s, gpu_s
    with JOBs_json.open() as fp:
        json_job_s = [path.Path(fn) for fn in json.load(fp)]
    for json_job in json_job_s:
        if not json_job.status(): continue
        job = Job.from_json(json_job)
        for method in job.task:
            task_s = job.get_task(method, status='RUN')
            for index, task in task_s:
                host = task['resource'][1]
                is_gpu = task['resource'][2]
                if is_gpu:
                    host_name = host.split("/")[0]
                    gpu_id = int(host.split("/")[1])
                    if host_name not in gpu_s:
                        gpu_s[host_name] = []
                    gpu_s[host_name].append(gpu_id)
                else:
                    cpu_s.append(host)
    return cpu_s, gpu_s

def assign_resource(job, updated):
    has_new_job = False
    for method in job.task:
        task_s = job.get_task(method, status='WAIT')
        if len(task_s) > 0:
            has_new_job = True
            break
    if not has_new_job:
        return updated

    if HOSTs_json.status():
        with HOSTs_json.open() as fp:
            host_s = json.load(fp)
    else:
        host_s = {}
    #
    gpu_host_s = {}
    for host in host_s:
        if not host_s[host][0]:
            continue
        gpu = GPU(host=host, login=HOSTNAME)
        gpu._visible_devices = host_s[host][1]
        gpu_host_s[host] = gpu
    #
    cpu_taken, gpu_taken = get_resource_taken()
    #
    cpu_s = []
    gpu_s = []
    for host in host_s:
        if host not in cpu_taken:
            cpu_s.append(host)
        if not host_s[host][0]:
            continue
        gpu = gpu_host_s[host]
        usable = gpu.usable(update=True)
        if host in gpu_taken:
            usable = [x for x in usable if x not in gpu_taken[host]]
        if len(usable) >= 0:
            for gpu_id in usable:
                gpu_s.append('%s/%d'%(host, gpu_id))
    cpu_status = [True for _ in cpu_s]
    gpu_status = [True for _ in gpu_s]
    #
    for method in job.task:
        task_s = job.get_task(method, status='WAIT')
        for index, task in task_s:
            if task['resource'][2]: # is GPU job
                if True not in gpu_status:
                    continue
                i = gpu_status.index(True)
                host = gpu_s[i]
                gpu_status[i] = False
            else:
                if True not in cpu_status:
                    continue
                i = cpu_status.index(True)
                host = cpu_s[i]
                cpu_status[i] = False
            #
            updated = True
            job.update_task_host(method, index, host)
            job.update_task_status(method, index, 'RUN')
    return updated

def submit_task(job):
    pass

def run(job, wait_after_run):
    while True:
        updated, status = check_status(job)
        if status:
            if updated: job.to_json()
            break
        #
        if RUNNER_METHOD == 'run':
            updated = assign_resource(job, updated)
        else:
            updated = submit_task(job, updated)
        #
        if updated: job.to_json()
        #
        if wait_after_run:
            time.sleep(30)
        else:
            break
    return status

def get_outputs(job, method):
    task_s = job.get_task(method, status='DONE')
    out_s = []
    for index, task in task_s:
        out_s.append(task['output'])
    return out_s

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

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    if arg.input_pdb is not None:
        arg.input_pdb = path.Path(arg.input_pdb)
    #
    # init
    job = import_module("init").prep(arg)

    # hybrid
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

    # locPREFMD
    import_module("locPREFMD").prep(job, job.init_pdb)
    if not run(job, arg.wait_after_run):
        return 
    locPREFMD_out = get_outputs(job, "locPREFMD")[:len(job.init_pdb)]

    # define topology
    import_module("define_topology").prep(job, locPREFMD_out[0][0])

    # equil
    import_module("equil").prep(job, 0, [out[0] for out in locPREFMD_out], path.Path("%s/equil.json"%(DEFAULT_HOME)))
    if not run(job, arg.wait_after_run):
        return 
    
    # prod
    for i in range(len(job.init_pdb)):
        import_module("prod").prep(job, i, i, path.Path("%s/prod.json"%DEFAULT_HOME), 5)
    if not run(job, arg.wait_after_run):
        return
    prod_out = get_outputs(job, 'prod')
    
    # score
    import_module("score").prep(job, [out[0] for out in prod_out])
    if not run(job, arg.wait_after_run):
        return

    # average
    if len(job.init_pdb) == 1:
        import_module("average").prep(job, '%s.prod_0'%job.title, [0], path.Path("%s/average.json"%DEFAULT_HOME), rule='score')
        import_module("average").prep(job, '%s.prod_0.cluster'%job.title, [0], path.Path("%s/average.json"%DEFAULT_HOME), rule='cluster')
        if not run(job, arg.wait_after_run):
            return
    else:
        import_module("average").prep(job, job.title, [i for i in range(len(job.init_pdb))], path.Path("%s/average.json"%DEFAULT_HOME), rule='casp12')
        for i in range(len(job.init_pdb)):
            import_module("average").prep(job, '%s.prod_%d'%(job.title, i), [i], path.Path("%s/average.json"%DEFAULT_HOME), rule='score')
        if not run(job, arg.wait_after_run):
            return
    average_out = get_outputs(job, 'average')[:N_MODEL]

    # model
    job.work_home.chdir()
    model_home = job.work_home.subdir("model", build=True)
    prep_s = []
    for i,out in enumerate(average_out):
        prep_fn = model_home.fn("prep_%d.pdb"%(i+1))
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
    locPREFMD_out = get_outputs(job, "locPREFMD")[-N_MODEL:]
    #
    model_s = []
    for i,out in enumerate(locPREFMD_out):
        model_fn = model_home.fn("model_%d.pdb"%(i+1))
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
        system(['cp', out[0].short(), pdb_fn.short()], verbose=arg.verbose)
    #
    job.remove_from_joblist()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit()

