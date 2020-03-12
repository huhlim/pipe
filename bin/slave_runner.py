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
    def __init__(self, gpu, time_interval, verbose):
        self.task_s = []
        self.gpu = gpu
        self.time_interval = time_interval
        self.verbose = verbose
    def wait(self):
        time.sleep(self.time_interval)
    def update_tasks(self):
        if not JOBs_json.status():
            return
        with JOBs_json.open() as fp:
            json_job_s = [path.Path(fn) for fn in json.load(fp)]
            #
        for json_job in json_job_s:
            if not json_job.status(): continue
            job = Job.from_json(json_job)
            for method in job.task:
                task = job.get_task(method, host=HOSTNAME, status='RUN')
                for index,X in task:
                    is_gpu_job = X['resource'][2]
                    t = ['WAIT', json_job, method, index, is_gpu_job]
                    is_new = True
                    for t0 in self.task_s:
                        if (t0[1] == t[1]) and (t0[2] == t[2]) and (t0[3] == t[3]):
                            is_new = False ; break
                    if is_new:
                        self.task_s.append(t)
    def check_resources(self, proc_s):
        cpu_status = True
        gpu_status = {} ; self.gpu.check()
        for gpu_id in self.gpu.usable(update=True):
            gpu_status[gpu_id] = True
        for proc in proc_s:
            _, use_gpu, gpu_id, _ = proc
            if use_gpu:
                gpu_status[gpu_id] = False
            else:
                cpu_status = False
        resource_available = cpu_status
        for gpu_id in gpu_status:
            if gpu_status[gpu_id]:
                resource_available = True
                break
        return resource_available, cpu_status, gpu_status
    def run(self):
        proc_s = []
        while True:
            resource_available, cpu_status, gpu_status = self.check_resources(proc_s)
            #
            if resource_available:
                self.update_tasks()
                #
                for task_id,task in enumerate(self.task_s):
                    if task[0] == 'RUNNING': continue
                    if task[0] == 'FINISHED': continue
                    #
                    use_gpu = task[4]
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
                        sys.stdout.write("PROC: %s\n"%cmd)
                    proc = sp.Popen(cmd, shell=True)
                    time.sleep(1.)
                    #
                    if proc.poll() is None:
                        task[0] = 'RUNNING'
                        proc_s.append((task_id, use_gpu, gpu_id, proc))
                    else:
                        task[0] = 'FINISHED'
            #
            proc_status, proc_running = self.check_proc_status(proc_s)
            for status, proc in zip(proc_status, proc_s):
                if not status: continue
                task_id, use_gpu, gpu_id, _ = proc
                self.task_s[task_id][0] = 'FINISHED'
                if use_gpu:
                    gpu_status[gpu_id] = True
                else:
                    cpu_status = True
            proc_s = proc_running
            #
            self.wait()

    def get_cmd(self, task, use_gpu, gpu_id=None):
        wrt = []
        if use_gpu and gpu_id is not None:
            wrt.append("export CUDA_VISIBLE_DEVICES=%s"%gpu_id)

        json_job = task[1]
        method = task[2]
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
            fout.write(json.dumps(host_s, indent=2))
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
        fout.write(json.dumps(host_s, indent=2))

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
    arg = arg.parse_args()
    #
    if arg.use_gpu:
        gpu = register_host(GPU())
    else:
        gpu = None
    #
    Queue(gpu, arg.time_interval, arg.verbose).run()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        unregister_host()
        sys.exit()
    finally:
        unregister_host()

