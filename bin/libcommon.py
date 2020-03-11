#!/usr/bin/env python

import os
import sys
import json
import path
from string import digits
import subprocess as sp

HOSTNAME=os.getenv('HOSTNAME')

WORK_HOME = os.getenv("PREFMD_HOME")
assert WORK_HOME is not None

BIN_HOME = '%s/bin'%WORK_HOME
EXEC_HOME = '%s/exec'%BIN_HOME
DEFAULT_HOME = '%s/default'%BIN_HOME
HOST_HOME = '%s/hosts'%BIN_HOME
JOBs_json = path.Path("%s/job_s.json"%HOST_HOME)
HOSTs_json = path.Path("%s/host_s.json"%HOST_HOME)

RUNNER_METHOD = 'run'
#RUNNER_METHOD = 'submit'

CHARMMEXEC = '/home/huhlim/apps/charmm/current/charmm'
if os.path.exists("/usr/bin/mpirun"):
    MPI_COMMAND = "/usr/bin/mpirun -np 4"
    CHARMMEXEC_MPI = '%s /home/huhlim/apps/charmm/current/charmm.mpi'%(MPI_COMMAND)
    CHARMM_MPI = True
else:
    CHARMMEXEC_MPI = CHARMMEXEC
    CHARMM_MPI = False

MAX_ERROR = 20

def system(cmd, verbose=True, stdout=False, stdin=None, outfile=None, errfile=None, redirect=False):
    if type(cmd) == type(""):
        cmd = cmd.strip().split()
    if verbose:
        sys.stdout.write("CMD: " + " ".join(cmd) + '\n')
    #
    if (not stdout) and (outfile is None) and (errfile is None):
        sp.call(cmd, stdin=stdin)
    else:
        if errfile == '/dev/null':
            errfile = sp.DEVNULL
        if redirect:
            errfile = sp.STDOUT
        try:
            out = sp.check_output(cmd, stdin=stdin, stderr=errfile)
            if sys.version_info.major == 3:
                out = out.decode("utf8")
            if outfile is not None:
                outfile.write(out)
        except sp.CalledProcessError:
            return ''
        return out

def JSONserialize(X):
    if isinstance(X, path.Dir):
        return 'Dir %s'%X.short()
    elif isinstance(X, path.Path):
        return 'Path %s'%X.short()
    else:
        return X

def JSONdeserialize(X):
    if isinstance(X, list):
        out = []
        for Xi in X:
            out.append(JSONdeserialize(Xi))
        return out
    elif isinstance(X, dict):
        out = {}
        for Xi in X:
            out[Xi] = JSONdeserialize(X[Xi])
        return out
    else:
        if isinstance(X, str):
            x = X.split()
            if x[0] == 'Dir':
                X = path.Dir(x[1])
            elif x[0] == 'Path':
                X = path.Path(x[1])
        return X

class Job(dict):
    def __init__(self, work_dir='.', title='', build=False):
        super().__init__()
        self.title = title
        if build:
            self.work_home = path.Dir("%s/%s"%(work_dir, title), build=True)
            self.json_job = self.work_home.fn("job.json")
        self.task = {}
    def __repr__(self):
        return self.work_home.path()
    def has(self, key):
        return hasattr(self, key)
    def __getitem__(self, key):
        return getattr(self, key)
    def __setitem__(self, key, item):
        setattr(self, key, item)
    @classmethod
    def from_json(cls, fn):
        with fn.open() as fp:
            X = json.load(fp)
        #
        cwd = os.getcwd()
        work_home = fn.dirname()
        work_home.chdir()
        #
        job = cls()
        for key in X:
            job[key] = JSONdeserialize(X[key])
        #
        os.chdir(cwd)
        return job

    def to_json(self):
        for method in self.task:
            if len(self.task[method]) == 0:
                del (self.task[method])
        #
        cwd = os.getcwd()
        self.work_home.chdir()
        with self.json_job.open('wt') as fout:
            fout.write(json.dumps(self.__dict__, indent=2, default=JSONserialize))
        os.chdir(cwd)

    def add_task(self, method, input_arg, output_arg, *arg, **kwarg):
        if method not in self.task:
            self.task[method] = []
        #
        def is_same(X, Y):
            if isinstance(X, list):
                status = [is_same(Xi, Yi) for Xi,Yi in zip(X, Y)]
                return (not False in status)
            else:
                return X == Y

        #
        prev_exists = False
        for prev in self.task[method]:
            input_status = True
            for i in range(len(prev['input'])):
                if not is_same(prev['input'][i], input_arg[i]):
                    input_status = False ; break
            if input_status:
                prev_exists = True
            output_status = True
            for i in range(len(prev['output'])):
                if not is_same(prev['output'][i], output_arg[i]):
                    output_status = False ; break
            if output_status:
                prev['resource'][0] = 'DONE'
        if prev_exists: return
        #
        # status allocated_resource use_gpu n_proc
        self.task[method].append({\
                 'resource': ["WAIT", None, kwarg.get("use_gpu", False), kwarg.get("n_proc", 1)],\
                 'input': input_arg, \
                 'output': output_arg, \
                 'etc': arg\
                 })
                 
    def get_task(self, method, host=None, status=None, not_status=None):
        if method not in self.task:
            return []
        #
        task_s = []
        for i,task in enumerate(self.task[method]):
            if status is not None and task['resource'][0] != status:
                continue
            if not_status is not None and task['resource'][0] == not_status:
                continue
            if host is not None and task['resource'][1] != host:
                continue
            task_s.append((i, task))
        return task_s
    def update_task_status(self, method, index, status):
        self.task[method][index]['resource'][0] = status
    def update_task_host(self, method, index, host):
        self.task[method][index]['resource'][1] = host
    def write_submit_script(self, method, index, cmd, submit=False):
        pass
    def append_to_joblist(self):
        if JOBs_json.status():
            with JOBs_json.open() as fp:
                job_s = json.load(fp)
        else:
            job_s = []
        #
        if self.json_job.path() not in job_s:
            job_s.append(self.json_job.path())
        with JOBs_json.open('wt') as fout:
            fout.write(json.dumps(job_s, indent=2))
    def remove_from_joblist(self):
        if not JOBs_json.status():
            return
        with JOBs_json.open() as fp:
            job_s = json.load(fp)
        if self.json_job.path() in job_s:
            job_s.remove(self.json_job.path())
        with JOBs_json.open('wt') as fout:
            fout.write(json.dumps(job_s, indent=2))

