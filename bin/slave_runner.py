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

class Task(object):
    def __init__(self, id, json_job, method, index, gpu, output):
        self.status = 'WAIT'
        self.id = id
        self.json_job = json_job
        self.method = method
        self.index = index
        self.gpu = gpu
        self.output = output
    def __repr__(self):
        wrt = []
        wrt.append("%d"%self.id)
        wrt.append(self.status)
        wrt.append(self.json_job.path())
        wrt.append(self.method)
        wrt.append("%d"%self.index)
        if self.gpu[0]:
            wrt.append("GPU/%d"%self.gpu[1])
        else:
            wrt.append("CPU  ")
        return ' '.join(wrt)
    def __eq__(self, othr):
        if self.json_job != othr.json_job:
            return False
        elif self.method != othr.method:
            return False
        elif self.index != othr.index:
            return False
        elif self.gpu[0] != othr.gpu[0]:
            return False
        elif self.gpu[0] and (self.gpu[1] != othr.gpu[1]):
            return False
        return True
    def check_output(self):
        for out in self.output:
            if not out.status():
                return False
        return True
    def update_status(self, status):
        self.status = status
    def to_json(self):
        return self.__repr__()
    @classmethod
    def from_json(cls, X, prev_s):
        X = X.split()
        if X[5].startswith("GPU/"):
            gpu = (True, int(X[5].split("/")[1]))
        else:
            gpu = (False, None)
        task = cls(int(X[0]), path.Path(X[2]), X[3], int(X[4]), gpu, None)
        task.update_status(X[1])
        if task in prev_s:
            prev = prev_s[prev_s.index(task)]
            task.output = prev.output
            return task
        #
        job = Job.from_json(task.json_job)
        task_s = job.get_task(task.method, host=HOSTNAME, status='RUN')
        for index,task in task_s:
            if index == task.index:
                task.output = task['output']
                break
        if task.output is None:
            sys.stderr.write("ERROR: failed to retrieve a task, %s\n"%task)
            return None
        return task

class Queue(object):
    def __init__(self, host_json, gpu, time_interval, verbose):
        self.host_json = host_json
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
                    if is_gpu_job:
                        gpu_id = int(X['resource'][1].split("/")[-1])
                    else:
                        gpu_id = None
                    t = Task(len(self.task_s), \
                            json_job, method, index, (is_gpu_job, gpu_id), X['output'])
                    if t not in self.task_s:
                        self.task_s.append(t)
        for i,task in enumerate(self.task_s):
            if task.status != 'CHECK': continue
            if task.check_output():
                task.update_status("FINISHED")
            else:
                task.update_status("WAIT")
    def from_json(self):
        task_s = []
        with self.host_json.open("r") as fp:
            for X in json.load(fp):
                task = Task.from_json(X, self.task_s)
                if task is not None:
                    task_s.append(task)
        return task_s
    def to_json(self):
        if self.host_json.status():
            comm_s = self.from_json()
            for task in self.task_s:
                if task.status in ['FINISHED']:
                    continue
                if task not in comm_s:
                    continue
                comm = comm_s[comm_s.index(task)]
                if comm.status in ['TERMINATE']:
                    task.update_status(comm.status)

        with self.host_json.open("wt") as fout:
            fout.write(json.dumps([task.to_json() for task in self.task_s], indent=2))

    def check_resources(self, proc_s):
        cpu_status = True
        if self.gpu is not None:
            gpu_status = {} ; self.gpu.check()
            for gpu_id in self.gpu.usable(update=True):
                gpu_status[gpu_id] = True
        else:
            gpu_status = {}
        for _,task in proc_s:
            use_gpu, gpu_id = task.gpu
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
                try:
                    self.update_tasks()
                except:
                    pass
                #
                for task in self.task_s:
                    if task.status in ['RUNNING', 'CHECK', 'FINISHED']:
                        continue
                    #
                    cmd = self.get_cmd(task)
                    use_gpu, gpu_id = task.gpu
                    if use_gpu:
                        if gpu_id not in gpu_status: continue
                        if not gpu_status[gpu_id]: continue
                        gpu_status[gpu_id] = False
                        #
                        env = os.environ.copy()
                        env['CUDA_VISIBLE_DEVICES'] = '%d'%gpu_id
                        proc = sp.Popen(cmd, env=env)
                    else:
                        if not cpu_status: continue
                        cpu_status = False
                        proc = sp.Popen(cmd)
                    if self.verbose:
                        sys.stdout.write("PROC: %s\n"%cmd)
                    time.sleep(1.)
                    #
                    if proc.poll() is None:
                        task.update_status("RUNNING")
                        proc_s.append((proc, task))
                    else:
                        task.update_status("CHECK")
            #
            proc_status, proc_running = self.check_proc_status(proc_s)
            for status, proc in zip(proc_status, proc_s):
                if not status: continue
                #
                _,task = proc
                task.update_status("CHECK")
                if task.gpu[0]:
                    gpu_status[task.gpu[1]] = True
                else:
                    cpu_status = True
            #
            self.to_json()
            #
            proc_s = []
            for proc,task in proc_running:
                t = self.task_s[self.task_s.index(task)]
                if t.status == 'TERMINATE':
                    print ("killing... %d"%proc.pid)
                    proc.terminate()
                    proc.wait()
                else:
                    proc_s.append((proc, task))
            for task in self.task_s:
                if task.status == 'TERMINATE':
                    task.update_status("FINISHED")
            #
            self.wait()

    def get_cmd(self, task):
        return ["%s/%s.py"%(BIN_HOME, task.method), 'run', task.json_job.short()]
    def check_proc_status(self, proc_s):
        status = [] ; proc_running = []
        for proc in proc_s:
            if proc[0].poll() is not None:
                status.append(True)
            else:
                status.append(False)
                proc_running.append(proc)
        return status, proc_running

def register_host(gpu):
    host_json = path.Path("%s/%s.json"%(HOST_HOME, HOSTNAME))
    #
    if HOSTs_json.status():
        with HOSTs_json.open() as fp:
            host_s = json.load(fp)
    else:
        host_s = {}
    #
    if HOSTNAME not in host_s:
        if gpu is None:
            host_s[HOSTNAME] = (False, [])
        else:
            host_s[HOSTNAME] = (True, gpu.CUDA_VISIBLE_DEVICES)
        with HOSTs_json.open('wt') as fout:
            fout.write(json.dumps(host_s, indent=2))
    else:
        sys.exit("%s is already in the host list\n"%HOSTNAME)
    return host_json, gpu

def unregister_host():
    host_json = path.Path("%s/%s.json"%(HOST_HOME, HOSTNAME))
    if host_json.status():
        host_json.remove()
    if not HOSTs_json.status():
        return
    with HOSTs_json.open() as fp:
        host_s = json.load(fp)
    if HOSTNAME in host_s:
        del host_s[HOSTNAME]
    with HOSTs_json.open("wt") as fout:
        fout.write(json.dumps(host_s, indent=2))

def main():
    arg = argparse.ArgumentParser(prog='slave_runner')
    arg.add_argument('--gpu', dest='use_gpu', action='store_true', default=False)
    arg.add_argument('--interval', dest='time_interval', default=60, type=int)
    arg.add_argument('--verbose', dest='verbose', action='store_true', default=False)
    #
    arg = arg.parse_args()
    #
    if arg.use_gpu:
        gpu = GPU()
        gpu.check()
        gpu.update()
        host_json, gpu = register_host(gpu)
    else:
        host_json, gpu = register_host(None)
    #
    Queue(host_json, gpu, arg.time_interval, arg.verbose).run()
    unregister_host()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        unregister_host()
        sys.exit()
    finally:
        unregister_host()
