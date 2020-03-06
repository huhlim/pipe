#!/usr/bin/env python

import os
import sys
import path
import json
import time
import argparse
import subprocess as sp

from libgpu import GPU
from libcommon import *

class Queue(object):
    def __init__(self, task_s, gpu, verbose):
        self.task_s = task_s
        self.gpu = gpu
        self.verbose = verbose
    def run(self):
        cpu_status = True
        gpu_status = {}
        for gpu_id in self.gpu.usable():
            gpu_status[gpu_id] = True
        #
        proc_s = []
        submitted = [False for _ in self.task_s]
        finished = [False for _ in self.task_s]
        while (False in finished):
            for task_id,task in enumerate(self.task_s):
                if finished[task_id]: continue
                if submitted[task_id]: continue
                #
                use_gpu = task[0]
                if use_gpu:
                    gpu_id_assigned = None
                    for gpu_id in gpu_status:
                        if gpu_status[gpu_id]:
                            gpu_id_assigned = gpu_id
                            break
                    if gpu_id_assigned is None:
                        continue
                    cmd = self.get_cmd(task, use_gpu, gpu_id=gpu_id)
                    gpu_status[gpu_id] = False
                else:
                    if not cpu_status: continue
                    cmd = self.get_cmd(task, use_gpu)
                    cpu_status = False
                    gpu_id = None
                #
                if self.verbose:
                    sys.stdout.write("%s\n"%cmd)
                proc = sp.Popen(cmd, shell=True)
                time.sleep(1.)
                if proc.poll() is None:
                    proc_s.append((task_id, use_gpu, gpu_id, proc))
                else:
                    finished[task_id] = True
                submitted[task_id] = True
            #
            proc_status, proc_running = self.check_proc_status(proc_s)
            for status, proc in zip(proc_status, proc_s):
                if not status: continue
                task_id, use_gpu, gpu_id, _ = proc
                finished[task_id] = True
                if use_gpu:
                    gpu_status[gpu_id] = True
                else:
                    cpu_status = True

            proc_s = proc_running
            time.sleep(30.)

    def get_cmd(self, task, use_gpu, gpu_id=None):
        wrt = []
        if use_gpu and gpu_id is not None:
            wrt.append("export CUDA_VISIBLE_DEVICES=%s"%gpu_id)
        use_gpu, method, json_job = task
        wrt.append("%s/%s.py run %s"%(BIN_HOME, method, json_job.short()))
        wrt = ' ; '.join(wrt)
        return wrt
    def check_proc_status(self, proc_s):
        status = [] ; proc_running = []
        for proc in proc_s:
            if proc[-1].poll() is not None:
                status.append(True)
            else:
                status.append(False)
                proc_running.append(proc)
        return status, proc_running

def register_host(gpu):
    if HOSTs_json.status():
        with HOSTs_json.open() as fp:
            host_s = json.load(fp)
    else:
        host_s = {}
    #
    if HOSTNAME not in host_s:
        host_s[HOSTNAME] = (True, gpu.CUDA_VISIBLE_DEVICES)
        with HOSTs_json.open('wt') as fout:
            fout.write(json.dumps(host_s))
    else:
        host = host_s[HOSTNAME]
        if host[0]:
            gpu._visible_devices = host[1]
        else:
            gpu = None
    return gpu

def unregister_host():
    if not HOSTs_json.status():
        return
    with HOSTs_json.open() as fp:
        host_s = json.load(fp)
    if HOSTNAME in host_s:
        del host_s[HOSTNAME]
    with HOSTs_json.open("wt") as fout:
        fout.write(json.dumps(host_s))

def get_tasks(json_job_s):
    task_s = []   # cpu / gpu
    #
    for json_job in json_job_s:
        if not json_job.status(): continue
        job = Job.from_json(json_job)
        for method in job.task:
            task = job.get_task(method, host=HOSTNAME, status='RUN')
            for _,X in task:
                is_gpu_job = X['resource'][2]
                t = (is_gpu_job, method, json_job)
                task_s.append((is_gpu_job, method, json_job))
    return task_s

def main():
    arg = argparse.ArgumentParser(prog='slave_runner')
    arg.add_argument('--gpu', dest='use_gpu', action='store_true', default=False)
    arg.add_argument('--interval', dest='time_interval', default=60, type=int)
    arg.add_argument('--verbose', dest='verbose', action='store_true', default=False)
    #
    if not JOBs_json.status():
        sys.stderr.write("There is no available job\n")
        return
    #
    arg = arg.parse_args()
    #
    if arg.use_gpu:
        gpu = register_host(GPU())
    else:
        gpu = None
    #
    while True:
        with JOBs_json.open() as fp:
            job_s = [path.Path(fn) for fn in json.load(fp)]
        if len(job_s) == 0: 
            time.sleep(arg.time_interval)
            continue
        #
        task_s = get_tasks(job_s)
        Queue(task_s, gpu, arg.verbose).run()
        #
        time.sleep(arg.time_interval)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        unregister_host()
        sys.exit()
    finally:
        unregister_host()

