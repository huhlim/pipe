#!/usr/bin/env python

import os
import sys
import json
import path
import signal
from string import digits
import subprocess as sp

assert sys.version_info.major == 3

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

SUBMIT_TEMPLATE = '%s/SUBMIT_TEMPLATE'%DEFAULT_HOME

CHARMMEXEC = '/home/huhlim/apps/charmm/current/charmm'
if os.path.exists("/usr/bin/mpirun"):
    MPI_COMMAND = "/usr/bin/mpirun -np 4"
    CHARMMEXEC_MPI = '%s /home/huhlim/apps/charmm/current/charmm.mpi'%(MPI_COMMAND)
    CHARMM_MPI = True
else:
    CHARMMEXEC_MPI = CHARMMEXEC
    CHARMM_MPI = False

N_MODEL = 5
MAX_ERROR = 20

MAX_SUBMIT = 10000
USERNAME = "heolim"

#TBM_EXCLUDE = '%s/exclude.casp13'%DEFAULT_HOME
TBM_EXCLUDE = None
HH_sequence_database = "/green/s2/huhlim/db/hhsuite/uc30/current/uc30"
HH_pdb70_database = "/green/s2/huhlim/work/database/hhsuite/pdb70/current/pdb70"

class GracefulExit(Exception):
    pass
def listen_signal():
    def gracefulExit(*arg):
        raise GracefulExit()
    signal.signal(signal.SIGTERM, gracefulExit)

def system(cmd, verbose=True, stdout=False, stdin=None, outfile=None, errfile=None):
    if type(cmd) == type(""):
        cmd = cmd.strip().split()
    if verbose:
        sys.stdout.write("CMD: " + " ".join(cmd) + '\n')
    #
    if stdout or (outfile is not None):
        STDOUT = sp.PIPE
    else:
        STDOUT = None
    if errfile is None:
        STDERR = None
    elif errfile == '/dev/null':
        STDERR = sp.DEVNULL
    else:
        STDERR = errfile
    #
    listen_signal()
    #
    proc_output = None
    try:
        proc = sp.Popen(cmd, stdin=stdin, stdout=STDOUT, stderr=STDERR)
        #print (proc.pid, ' '.join(cmd))
        proc_output = proc.communicate()
    except (GracefulExit, KeyboardInterrupt) as error:
        #print ("killing... %d"%proc.pid)
        proc.terminate()
        sys.exit()
    #
    if (proc_output is not None) and (STDOUT is not None):
        out = proc_output[0].decode("utf8")
        if outfile is not None:
            outfile.write(out)
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
            if RUNNER_METHOD == 'submit':
                self.queue_home = self.work_home.subdir("queue", build=True)
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
        if RUNNER_METHOD == 'submit':
            job.queue_home = job.work_home.subdir("queue", build=True)
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
            if host is not None and (task['resource'][1] != host and task['resource'][1].split("/")[0] != host):
                continue
            task_s.append((i, task))
        return task_s
    def update_task_status(self, method, index, status):
        self.task[method][index]['resource'][0] = status
    def update_task_host(self, method, index, host):
        self.task[method][index]['resource'][1] = host
    def write_submit_script(self, method, index, cmd_s):
        with open(SUBMIT_TEMPLATE) as fp:
            que = fp.read()
        job_name = '%s.%s.%d'%(self.title, method, index)
        que = que.replace("JOB_NAME", job_name)
        que = que.replace("JOB_OUTPUT", self.queue_home.fn("%s.log"%(job_name)).path())
        #
        que_fn = self.queue_home.fn("%s.sh"%job_name)
        with que_fn.open("wt") as fout:
            fout.write(que)
            fout.writelines(cmd_s)
        #
        sbatch = system(['sbatch', que_fn.short()], stdout=True)
        que_id = sbatch.strip().split()[-1]

        self.update_task_status(method, index, "SUBMITTED")
        self.update_task_host(method, index, que_id)
        #
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

