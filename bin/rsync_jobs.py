#!/usr/bin/env python

import os
import sys
import path
import subprocess as sp

from libcommon import *

EXPAND_s = {}
EXPAND_s["hybrid"] = "model_s"
EXPAND_s["average"] = "pdb_s"


class FakeJob(object):
    def __init__(self, fn):
        self.work_home = path.Dir(fn)
        self.title = fn.split("/")[-1]


### THIS IS DIFFERENT FROM libmain.get_outputs
def get_outputs(job, method):
    task_s = job.get_task(method, status="DONE")
    if method in EXPAND_s:
        expand = EXPAND_s[method]
    else:
        expand = None
    #
    out_s = []
    for index, task in task_s:
        if expand is None:
            out_s.append(task["output"])
        else:
            output_expanded = []
            for out in task["output"]:
                if out.endswith(expand):
                    output_expanded.append(out)
                    with out.open() as fp:
                        _out = []
                        for line in fp:
                            if line.startswith("#"):
                                continue
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
    # output_s.append(job.json_job)
    #
    if job.run_type == "refine":
        output_s.extend(job.init_home.glob("*"))
    elif job.run_type == "refine_remote":
        output_s.extend(job.init_home.glob("*"))
    elif job.run_type == "sp":
        output_s.append(job.init_fa)
    #
    for method in job.task:
        out_s = get_outputs(job, method)
        output_s.extend(flatten_list(out_s))
    #
    if job.run_type == "sp" and job.has("refine_s"):
        for refine_home in job.refine_s:
            refine_json_job = refine_home.fn("job.json")
            refine_job = Job.from_json(refine_json_job)
            output_s.extend(get_job_file_list(refine_job))
    #
    if job.run_type in ["refine", "refine_remote"]:
        if job.has("is_oligomer") and job.work_home.subdir("oligomer").status():
            output_s.append(job.work_home.subdir("oligomer"))
        if job.has("is_membrane_protein") and job.work_home.subdir("membrane").status():
            output_s.append(job.work_home.subdir("membrane"))
        if job.has("has_ligand") and job.work_home.subdir("ligand").status():
            output_s.append(job.work_home.subdir("ligand"))
    #
    output_s.extend(job.work_home.subdir("model").glob("*"))
    output_s.extend(job.work_home.subdir("final").glob("*"))
    return output_s


def receive_via_rsync(remote, job):
    with open("SEND", "wt") as fout:
        hostname, work_home = remote.split(":", maxsplit=1)
        cmd = [
            "ssh",
            hostname,
            "$PIPE_HOME/bin/rsync_jobs.py output %s/%s/job.json" % (work_home, job.title),
        ]
        output_s = sp.check_output(cmd).decode("utf8")
        fout.write(output_s)

    cmd = []
    cmd.append("rsync")
    cmd.append("-ard")
    cmd.append("--files-from=SEND")
    cmd.append("%s/%s/" % (remote, job.title))
    cmd.append("./")
    sp.call(cmd)


def send_via_rsync(remote, job, output_s):
    flist = open("SEND", "wt")
    for fn in output_s:
        flist.write("%s\n" % fn.short())
    flist.close()
    #
    cmd = []
    cmd.append("rsync")
    cmd.append("-ard")
    cmd.append("--files-from=SEND")
    cmd.append("./")
    cmd.append("%s/%s/" % (remote, job.title))
    sp.call(cmd)


def main():
    if len(sys.argv) < 3:
        sys.stderr.write("usage: %s [METHOD] [REMOTE] [JOB]\n" % __file__)
        return
    #
    method = sys.argv[1]
    if method not in ["send", "receive", "output"]:
        sys.stderr.write("error: [METHOD = (send/receive)]\n")
        return
    if method != "output":
        remote = sys.argv[2]
        json_job = path.Path(sys.argv[3])
    else:
        json_job = path.Path(sys.argv[2])
    #
    if json_job.endswith("job.json"):
        job = Job.from_json(json_job)
    elif method != "receive":
        sys.exit("job.json MUST BE PROVIDED")
    else:
        job = FakeJob(json_job)
    #
    cwd = os.getcwd()
    job.work_home.chdir()
    #
    if method == "output":
        output_s = get_job_file_list(job)
        for output in output_s:
            sys.stdout.write("%s\n" % output.short())
    elif method == "send":
        output_s = get_job_file_list(job)
        send_via_rsync(remote, job, output_s)
    else:
        receive_via_rsync(remote, job)


if __name__ == "__main__":
    main()
