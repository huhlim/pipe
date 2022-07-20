#!/usr/bin/env python

import os
import sys
import path
import json
import argparse
import mdtraj

from libcommon import *

METHOD = "prod_meta"
EXEC = "%s/prod_meta.py" % EXEC_HOME


def prep(job, prod_index, prev_job):
    # if len(job.get_task(METHOD, not_status='DONE')) > 0:
    #    return
    #
    job.prod_home = job.work_home.subdir("prod", build=True)
    job.prod_home.chdir()
    #
    iter_home = job.prod_home.subdir("%d" % prod_index, build=True)
    top_index_fn = job.top_index_s[prod_index]
    #
    prev_tasks = prev_job.get_task("prod")
    for i, prev in prev_tasks:
        run_home = iter_home.subdir("%d" % i)
        #
        input_s = [run_home, top_index_fn, prev_job.top_fn, prev["output"][0]]
        output_s = [run_home.fn("solute.dcd")]
        job.add_task(METHOD, input_s, output_s, n_proc=12)
    #
    job.to_json()


def run(job):
    task_s = job.get_task(METHOD, host=HOSTNAME, status="RUN")
    if len(task_s) == 0:
        return
    #
    for index, task in task_s:
        run_home = task["input"][0]
        index_fn = task["input"][1]
        top_fn = task["input"][2]
        in_dcd = task["input"][3]
        #
        output_s = task["output"]
        status = True
        for output in output_s:
            if not output.status():
                status = False
                break
        if status:
            continue
        #
        run_home.build()
        run_home.chdir()
        #
        cmd = [EXEC]
        cmd.extend(["--index", index_fn.short()])
        cmd.extend(["--top", top_fn.short()])
        cmd.extend(["--dcd", in_dcd.short()])
        cmd.extend(["--out", output_s[0].short()])
        #
        system(cmd, verbose=job.verbose)


def submit(job):
    task_s = job.get_task(METHOD, status="SUBMIT")
    if len(task_s) == 0:
        return
    #
    for index, task in task_s:
        run_home = task["input"][0]
        index_fn = task["input"][1]
        top_fn = task["input"][2]
        in_dcd = task["input"][3]
        #
        output_s = task["output"]
        status = True
        for output in output_s:
            if not output.status():
                status = False
                break
        if status:
            continue
        #
        run_home.build()
        run_home.chdir()
        #
        cmd = [EXEC]
        cmd.extend(["--index", index_fn.short()])
        cmd.extend(["--top", top_fn.short()])
        cmd.extend(["--dcd", in_dcd.short()])
        cmd.extend(["--out", output_s[0].short()])
        cmd_s.append(" ".join(cmd) + "\n")
        #
        job.write_submit_script(METHOD, index, cmd_s)


def status(job):
    pass


def main():
    arg = argparse.ArgumentParser(prog="prod")
    arg.add_argument(dest="command", choices=["prep", "run"], help="exec type")
    arg.add_argument(dest="work_dir", help="work_dir, which has a JSON file")
    arg.add_argument("--index", dest="prod_index", help="prod_index", type=int)
    arg.add_argument(
        "-i",
        "--input",
        dest="input_equil",
        nargs="*",
        help='input equil index, mandatory for "prep"',
    )
    arg.add_argument(
        "-j", "--json", dest="input_json", help='input JSON file, mandatory for "prep"'
    )
    arg.add_argument("-n", dest="n_replica", type=int, default=5, help="number of replicas")

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
        if arg.input_json is None:
            sys.exit("Error: input_json required\n")
        arg.input_pdb = [path.Path(fn) for fn in arg.input_pdb]
        #
        prep(job, arg.prod_index, arg.input_equil, path.Path(arg.input_json), arg.n_replica)

    elif arg.command == "run":
        run(job)


if __name__ == "__main__":
    main()
