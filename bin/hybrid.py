#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

METHOD = "hybrid"
EXEC = "%s/hybrid.py" % EXEC_HOME


def prep(job, input_pdb):
    if len(job.get_task(METHOD, not_status="DONE")) > 0:
        return
    #
    job.hybrid_home = job.work_home.subdir("hybrid", build=True)
    out = job.hybrid_home.fn("model_s")
    #
    job.add_task(METHOD, [job.title, input_pdb, job.hybrid_home], [out], use_gpu=False, n_proc=40)
    #
    job.to_json()


def run(job):
    task_s = job.get_task(METHOD, host=HOSTNAME, status="RUN")
    if len(task_s) == 0:
        return
    #
    for index, task in task_s:
        title = task["input"][0]
        input_pdb = task["input"][1]
        run_home = task["input"][2]
        output_list = task["output"][0]
        if output_list.status():
            continue
        n_proc = task["resource"][3]
        #
        run_home.chdir()
        cmd = [EXEC, title, input_pdb.short(), "-j", "%d" % n_proc]
        system(cmd, verbose=job.verbose)


def submit(job):
    task_s = job.get_task(METHOD, status="SUBMIT")
    if len(task_s) == 0:
        return
    #
    for index, task in task_s:
        title = task["input"][0]
        input_pdb = task["input"][1]
        run_home = task["input"][2]
        output_list = task["output"][0]
        if output_list.status():
            continue
        n_proc = task["resource"][3]
        #
        run_home.chdir()
        #
        cmd = []
        cmd.append("cd %s\n" % run_home)
        cmd.append(" ".join([EXEC, title, input_pdb.short(), "-j", "%d" % n_proc]) + "\n")
        #
        job.write_submit_script(METHOD, index, cmd)


def status(job):
    task_s = job.get_task(METHOD)
    if len(task_s) == 0:
        return
    #
    for index, task in task_s:
        output_pdb = task["output"][0]
        if output_pdb.status():
            job.update_task_status(METHOD, index, "DONE")
        elif output_pdb.exists():
            job.update_task_status(METHOD, index, "RUN")


def main():
    arg = argparse.ArgumentParser(prog="hybrid")
    arg.add_argument(dest="command", choices=["prep", "run"], help="exec type")
    arg.add_argument(dest="work_dir", help="work_dir, which has a JSON file")
    arg.add_argument("-i", "--input", dest="input_pdb", help='input PDB file, mandatory for "prep"')

    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args()

    #
    if arg.work_dir.endswith(".json"):
        arg.json_job = path.Path(arg.work_dir)
    else:
        arg.work_dir = path.Dir(arg.work_dir)
        arg.json_job = arg.work_dir.fn("job.json")
    #
    job = Job.from_json(arg.json_job)
    #
    if arg.command == "prep":
        if arg.input_pdb is None:
            sys.exit("Error: input_pdb required\n")
        arg.input_pdb = path.Path(arg.input_pdb)
        #
        prep(job, arg.input_pdb)

    elif arg.command == "run":
        run(job)


if __name__ == "__main__":
    main()
