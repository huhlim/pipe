#!/usr/bin/env python

import os
import sys
import path
import subprocess as sp
from tempfile import TemporaryFile

from libcommon import *

ENDPOINT = {}
ENDPOINT['huhlim']  = 'b66c90e4-1379-11e8-b5f8-0ac6873fc732'
ENDPOINT['comet']   = 'de463f97-6d04-11e5-ba46-22000b92c6ec'
ENDPOINT['oasis']   = 'c982d486-6d04-11e5-ba46-22000b92c6ec'
ENDPOINT['hpcc']    = 'a640bafc-6d04-11e5-ba46-22000b92c6ec'

ENDPATH = {}
ENDPATH['huhlim']   = '~/'
ENDPATH['comet']    = '~/'
ENDPATH['oasis']    = '/oasis/scratch-comet/huhlim/temp_project'
ENDPATH['hpcc']     = '~/'

EXPAND_s = {}
EXPAND_s['hybrid'] = 'model_s'
EXPAND_s['average'] = 'pdb_s'

def check_output(cmd):
    return sp.check_output(cmd).decode("utf-8")

def initialize_globus():
    try:
        status = check_output(["globusconnect", '-status'])
    except sp.CalledProcessError as error:
        status = error.output.decode("utf-8")
    status = status.split("\n")[0]
    #
    if not status.startswith("No Globus "):
        return None
    #
    proc_globus = sp.Popen(["globusconnect", '-start'])
    return proc_globus

def finalize_globus(proc_globus):
    if proc_globus is None:
        return
    #
    try:
        status = check_output(["globusconnect", '-status'])
    except sp.CalledProcessError as error:
        status = error.output.decode("utf-8")
    status = status.split("\n")[0]
    #
    if not status.startswith("No Globus "):
        return
    #
    check_output(['globusconnect', '-stop'])

### THIS IS DIFFERENT FROM libmain.get_outputs
def get_outputs(job, method):
    task_s = job.get_task(method, status='DONE')
    if method in EXPAND_s:
        expand = EXPAND_s[method]
    else:
        expand = None
    #
    out_s = []
    for index, task in task_s:
        if expand is None:
            out_s.append(task['output'])
        else:
            output_expanded = []
            for out in task['output']:
                if out.endswith(expand):
                    output_expanded.append(out)
                    with out.open() as fp:
                        _out = []
                        for line in fp:
                            if line.startswith("#"): continue
                            _out.append(path.Path(line.strip()))
                        output_expanded.append(_out)
                else:
                    output_expanded.append(out)
            out_s.append(output_expanded)
    return out_s

def flatten_list(input_list):
    output_list = []
    for item in input_list:
        if isinstance(item, list):
            output_list.extend(flatten_list(item))
        else:
            output_list.append(item)
    return output_list

def get_job_file_list(job):
    output_s = []
    #output_s.append(job.json_job)
    #
    if job.run_type == 'refine':
        output_s.extend(job.init_home.glob("*"))
    elif job.run_type == 'sp':
        output_s.append(job.init_fa)
    #
    for method in job.task:
        out_s = get_outputs(job, method)
        output_s.extend(flatten_list(out_s))
    #
    if job.run_type == 'sp' and job.has("refine_s"):
        for refine_home in job.refine_s:
            refine_json_job = refine_home.fn("job.json")
            refine_job = Job.from_json(refine_json_job)
            output_s.extend(get_job_file_list(refine_job))
    return output_s

def set_remote_path(remote, title):
    host, remote_dir = remote.split(":")
    if host in ['oasis'] and remote_dir.startswith("~/"):
        remote_dir = '%s/%s'%(ENDPATH[host], remote_dir[2:])
    remote_path = '%s:%s/%s'%(ENDPOINT[host], remote_dir, title)
    return remote_path

def get_remote_file_list(job, remote):
    remote_path = set_remote_path(remote, job.title)
    #
    cmd = ['globus', 'ls', remote_path]
    output = check_output(cmd).split("\n")[:-1]
    #
    out_s = []
    for out in output:
        if out not in ['job.json', 'queue/']:
            out_s.append(out)
    return out_s

def send_via_globus(remote, job, out_s):
    remote_path = set_remote_path(remote, job.title)
    curr_path = '%s:%s'%(ENDPOINT['huhlim'], os.getcwd())
    #
    flist = TemporaryFile(mode='w+t')
    for out in out_s:
        fn = out.short()
        flist.write("%s %s\n"%(fn, fn))
    #
    cmd = ['globus', 'transfer', '--sync-level', 'mtime']
    cmd.append(curr_path)
    cmd.append(remote_path)
    cmd.append("--batch")
    #
    flist.seek(0)
    sys.stdout.write(" ".join(cmd) + '\n')
    sp.call(cmd, stdin=flist)
    #
    flist.close()

def receive_via_globus(remote, job, out_s):
    remote_path = set_remote_path(remote, job.title)
    curr_path = '%s:%s'%(ENDPOINT['huhlim'], os.getcwd())
    #
    flist = TemporaryFile(mode='w+t')
    for fn in out_s:
        if fn.endswith("/"):
            flist.write("%s %s --recursive\n"%(fn, fn))
        else:
            flist.write("%s %s\n"%(fn, fn))
    #
    cmd = ['globus', 'transfer', '--sync-level', 'mtime']
    cmd.append(remote_path)
    cmd.append(curr_path)
    cmd.append("--batch")
    #
    flist.seek(0)
    sys.stdout.write(" ".join(cmd) + '\n')
    sp.call(cmd, stdin=flist)
    #
    flist.close()

def main():
    if len(sys.argv) < 4:
        sys.stderr.write("usage: %s [METHOD] [REMOTE] [JOB]\n"%__file__)
        return
    #
    method = sys.argv[1]
    if method not in ['send', 'receive']:
        sys.stderr.write("error: [METHOD = (send/receive)]\n")
        return
    remote = sys.argv[2]
    json_job = path.Path(sys.argv[3])
    #
    job = Job.from_json(json_job)
    cwd = os.getcwd()
    job.work_home.chdir()
    #
    proc_globus = initialize_globus()
    #
    if method == 'send':
        output_s = get_job_file_list(job)
        send_via_globus(remote, job, output_s)
    else:
        output_s = get_remote_file_list(job, remote)
        receive_via_globus(remote, job, output_s)
        #
    finalize_globus(proc_globus)

if __name__ == '__main__':
    main()
