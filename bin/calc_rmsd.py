#!/usr/bin/env python

import os
import sys
import path
import json
import argparse

from libcommon import *

METHOD = "calc_rmsd"
EXEC = "%s/calc_tmscore" % EXEC_HOME


def prep(job, index, ref_fn_s, input_dcd):
    # if len(job.get_task(METHOD, not_status='DONE')) > 0:
    #    return
    #
    job.rmsd_home = job.work_home.subdir("rmsd", build=True)
    job.rmsd_home.chdir()
    #
    sub_home = job.rmsd_home.subdir("%d" % index, build=True)
    sub_home.chdir()
    #
    for i, dcd_fn in enumerate(input_dcd):
        if not dcd_fn.status():
            continue
        run_home = dcd_fn.dirname()
        input_s = [dcd_fn]
        output_s = [sub_home.fn("rmsd.%d.dat" % i)]
        job.add_task(METHOD, input_s, output_s, use_gpu=False, n_proc=24)
    #
    job.to_json()


def run(job):
    task_s = job.get_task(METHOD, host=HOSTNAME, status="RUN")
    if len(task_s) == 0:
        return
    #
    for index, task in task_s:
        input_dcd = task["input"][0]
        run_home = input_dcd.dirname()
        #
        output_qual = task["output"][0]
        if output_qual.status():
            continue
        #
        run_home.chdir()
        #
        pdblist = run_home.fn("pdb_s")

        if not output_qual.status():
            cmd = [EXEC]
            cmd.extend(["-j", "%d" % task["resource"][3]])
            cmd.extend(["-r", job.init_pdb[0].short()])
            cmd.extend(["-l", pdblist.short()])
            with output_qual.open("wt") as fout:
                system(cmd, outfile=fout, errfile="/dev/null")
    #


def submit(job):
    task_s = job.get_task(METHOD, status="SUBMIT")
    if len(task_s) == 0:
        return
    #
    for index, task in task_s:
        input_dcd = task["input"][0]
        run_home = input_dcd.dirname()
        #
        output_qual = task["output"][0]
        if output_qual.status():
            continue
        #
        run_home.chdir()
        #
        pdblist = run_home.fn("pdb_s")
        #
        cmd_s = []
        if not output_qual.status():
            cmd = [EXEC]
            cmd.extend(["-j", "%d" % task["resource"][3]])
            cmd.extend(["-r", job.init_pdb[0].short()])
            cmd.extend(["-l", pdblist.short()])
            cmd.append("> %s 2> /dev/null" % output_qual.short())
            cmd_s.append(" ".join(cmd) + "\n")
        #
        job.write_submit_script(METHOD, index, cmd_s)


def main():
    arg = argparse.ArgumentParser(prog="calc_rmsd")
    arg.add_argument(dest="command", choices=["prep", "run"], help="exec type")
    arg.add_argument(dest="work_dir", help="work_dir, which has a JSON file")
    arg.add_argument(
        "-i", "--input", dest="input_dcd", nargs="*", help='input DCD file, mandatory for "prep"'
    )

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
        if arg.input_dcd is None:
            sys.exit("Error: input_dcd required\n")
        arg.input_dcd = [path.Path(fn) for fn in arg.input_dcd]
        #
        prep(job, 0, arg.input_dcd)

    elif arg.command == "run":
        run(job)


if __name__ == "__main__":
    main()
