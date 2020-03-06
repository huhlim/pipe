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

CHARMMEXEC = '/home/huhlim/apps/charmm/current/charmm'
if os.path.exists("/usr/bin/mpirun"):
    MPI_COMMAND = "/usr/bin/mpirun -np 4"
    CHARMMEXEC_MPI = '%s /home/huhlim/apps/charmm/current/charmm.mpi'%(MPI_COMMAND)
    CHARMM_MPI = True
else:
    CHARMMEXEC_MPI = CHARMMEXEC
    CHARMM_MPI = False

def system(cmd, verbose=True, stdout=False, stdin=None, outfile=None, errfile=None, redirect=False):
    if type(cmd) == type(""):
        cmd = cmd.strip().split()
    if verbose:
        sys.stdout.write("CMD: " + " ".join(cmd) + '\n')
    #
    if (not stdout) and (outfile is None) and (errfile is None):
        sp.call(cmd, stdin=stdin)
    else:
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

class Job(dict):
    def __init__(self, work_dir='.', title='', build=False):
        super().__init__()
        self.title = title
        if build:
            self.work_home = path.Dir("%s/%s"%(work_dir, title), build=True)
            self.json_job = self.work_home.fn("job.json")
            #self.json_resource = self.work_home.fn("resource.json")
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
        def deserialize(X):
            if isinstance(X, list):
                out = []
                for Xi in X:
                    out.append(deserialize(Xi))
                return out
            elif isinstance(X, dict):
                out = {}
                for Xi in X:
                    out[Xi] = deserialize(X[Xi])
                return out
            else:
                if isinstance(X, str):
                    x = X.split()
                    if x[0] == 'Dir':
                        X = path.Dir(x[1])
                    elif x[0] == 'Path':
                        X = path.Path(x[1])
                return X

        with fn.open() as fp:
            X = json.load(fp)
        #
        cwd = os.getcwd()
        work_home = fn.dirname()
        work_home.chdir()
        #
        job = cls()
        for key in X:
            job[key] = deserialize(X[key])
        #
        os.chdir(cwd)
        return job

    def to_json(self):
        def serialize(X):
            if isinstance(X, path.Dir):
                return 'Dir %s'%X.short()
            elif isinstance(X, path.Path):
                return 'Path %s'%X.short()
            else:
                return X
        #
        for method in self.task:
            if len(self.task[method]) == 0:
                del (self.task[method])
        #
        cwd = os.getcwd()
        self.work_home.chdir()
        with self.json_job.open('wt') as fout:
            fout.write(json.dumps(self.__dict__, indent=2, default=serialize))
        os.chdir(cwd)

    def add_task(self, method, input_arg, output_arg, *arg, **kwarg):
        if method not in self.task:
            self.task[method] = []
        # status allocated_resource use_gpu n_proc
        self.task[method].append({\
                 'resource': ['WAIT', None, kwarg.get("use_gpu", False), kwarg.get("n_proc", 1)],\
                 'input': input_arg, \
                 'output': output_arg, \
                 'etc': arg\
                 })
                 
    def get_task(self, method, host=None, status=None):
        if method not in self.task:
            return []
        #
        task_s = []
        for i,task in enumerate(self.task[method]):
            if status is not None and task['resource'][0] != status:
                continue
            if host is not None and task['resource'][1].split("/")[0] != host:
                continue
            task_s.append((i, task))
        return task_s
    def update_task_status(self, method, index, status):
        self.task[method][index][0] = status
    def write_submit_script(self, method, index, cmd, submit=False):
        pass
