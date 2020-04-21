#!/usr/bin/env python

import path
import json
import time

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

def run(job, wait_after_run, sleep=30):
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
            time.sleep(sleep)
        else:
            break
    return status

def get_outputs(job, method):
    task_s = job.get_task(method, status='DONE')
    out_s = []
    for index, task in task_s:
        out_s.append(task['output'])
    return out_s
