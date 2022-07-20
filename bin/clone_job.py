#!/usr/bin/env python

import os
import sys
import path
import subprocess as sp

from libcommon import *

EXPAND_s = {}
EXPAND_s["hybrid"] = "model_s"
EXPAND_s["average"] = "pdb_s"


def check_output(cmd):
    return sp.check_output(cmd).decode("utf-8")


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
    return output_s


def main():
    if len(sys.argv) < 3:
        sys.stderr.write("usage: %s [JOB] [CLONE]\n" % __file__)
        sys.stderr.write("example: %s xxxx/job.json yyyy\n" % __file__)
        return
    #
    json_job = path.Path(sys.argv[1])
    clone_dir = path.Dir(sys.argv[2], build=True)
    clone_title = clone_dir.name()
    #
    job = Job.from_json(json_job)
    cwd = os.getcwd()
    job.work_home.chdir()
    #
    output_s = get_job_file_list(job)
    for output in output_s:
        fn = output.short()
        if job.title != clone_title:
            fn = fn.replace(job.title, clone_title)
        #
        clone_fn = clone_dir.fn(fn)
        clone_fn_dir = clone_fn.dirname()
        clone_fn_dir.build()
        #
        cmd = ["cp", "-v"]
        cmd.append(output.path())
        cmd.append(clone_fn.path())
        sp.call(cmd)


if __name__ == "__main__":
    main()
