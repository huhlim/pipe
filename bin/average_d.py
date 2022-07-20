#!/usr/bin/env python

import os
import sys
import path
import json
import argparse
import numpy as np

from libcommon import *

METHOD = "average_d"
EXEC = "%s/average_d.py" % EXEC_HOME

PARAM = {}
PARAM["score"] = ("RWplus", 25.0)
PARAM["casp12"] = ("RWplus", 0.5, 225.0, 45.0)
PARAM["cluster"] = (2.0, 20, 5)  # rmsd_cutoff, subsample, max_run_md
PARAM["msm"] = ("RWplus", 25.0)

PARAM_big = {}
PARAM_big["score"] = ("dDFIRE", 25.0)
PARAM_big["casp12"] = ("dDFIRE", 0.5, 225.0, 45.0)
PARAM_big["cluster"] = (2.0, 20, 5)  # rmsd_cutoff, subsample, max_run_md
PARAM_big["msm"] = ("dDFIRE", 25.0)


def prep(job, output_prefix, input_prod, input_json, rule="score"):
    #    if len(job.get_task(METHOD, not_status='DONE')) > 0:
    #        return
    #
    prod_s = job.get_task("prod")
    prod_s = [prod for _, prod in prod_s if int(prod["input"][0].split("/")[-2]) in input_prod]
    #
    score_s = job.get_task("score")
    #
    job.average_d_home = job.work_home.subdir("average_d", build=True)
    job.average_d_home.chdir()
    #
    if rule == "score":
        input_s = [output_prefix, (rule, PARAM[rule]), input_json, [], []]
    elif rule == "casp12":
        input_s = [output_prefix, (rule, PARAM[rule]), input_json, [], [], []]
    elif rule == "cluster":
        input_s = [output_prefix, (rule, PARAM[rule]), input_json, []]
    if job.n_atom > 22400:
        input_s[1] = (rule, PARAM_big[rule])
    #
    for prod in prod_s:
        if prod["resource"][0] != "DONE":
            return
        if not prod["output"][0].status():
            return

        if rule in ["score", "casp12"]:
            score = None
            for _, s in score_s:
                if s["input"][0] == prod["output"][0]:
                    score = s
                    break
            if score is None:
                return
            if score["resource"][0] != "DONE":
                return
            if not score["output"][0].status():
                return
        #
        input_s[3].append(prod["output"][0])
        if rule in ["score", "casp12"]:
            input_s[4].append(score["output"][0])
        if rule in ["casp12"]:
            input_s[5].append(score["output"][1])
    #
    output_s = [job.average_d_home.fn("%s.pdb_s" % output_prefix)]
    job.add_task(METHOD, input_s, output_s, use_gpu=True, n_proc=1)
    #
    job.to_json()


def prep_from_msm(job, output_prefix, msm_fn, input_json):
    job.average_d_home = job.work_home.subdir("average_d", build=True)
    job.average_d_home.chdir()
    #
    rule = "msm"
    input_s = [output_prefix, (rule, PARAM[rule]), input_json, [], [], msm_fn]
    #
    msm = np.load(msm_fn.short())
    traj_fn_s = [path.Path(fn) for fn in msm["traj_s"]]
    #
    score_s = job.get_task("score")
    for traj_fn in traj_fn_s:
        score = None
        for _, s in score_s:
            if s["input"][0] == traj_fn:
                score = s
                break
        if score is None:
            return
        if score["resource"][0] != "DONE":
            return
        if not score["output"][0].status():
            return
        #
        input_s[3].append(traj_fn)
        input_s[4].append(score["output"][0])
    #
    output_s = [job.average_d_home.fn("%s.pdb_s" % output_prefix)]
    job.add_task(METHOD, input_s, output_s, use_gpu=True, n_proc=1)
    #
    job.to_json()


def run(job):
    task_s = job.get_task(METHOD, host=HOSTNAME, status="RUN")
    if len(task_s) == 0:
        return
    gpu_id = os.environ["CUDA_VISIBLE_DEVICES"]
    #
    for index, task in task_s:
        if task["resource"][1].split("/")[1] != gpu_id:
            continue
        input_s = task["input"]
        output_prefix = input_s[0]
        rule = input_s[1]
        input_json = input_s[2]
        input_dcd_s = input_s[3]
        output_pdb = task["output"][0]
        if output_pdb.status():
            continue
        #
        with input_json.open() as fp:
            options = json.load(fp)
        options["ssbond"] = []
        for line in job.ssbond:
            chain_1 = line[15]
            chain_2 = line[29]
            if chain_1 == " " and chain_2 == " ":
                line = "%sA%sA%s" % (line[:15], line[16:29], line[30:])
            options["ssbond"].append(line)
        options["rule"] = rule
        #
        job.average_d_home.chdir()
        #
        input_json = job.average_d_home.fn("%s.json" % output_prefix)
        with input_json.open("wt") as fout:
            fout.write(json.dumps(options, indent=2, default=JSONserialize))
        #
        cmd = [EXEC, output_prefix, job.top_fn.short()]
        cmd.extend(["--input", input_json.short()])
        if rule[0] == "msm":
            cmd.extend(["--msm", input_s[5].short()])
        cmd.append("--dcd")
        cmd.extend([fn.short() for fn in input_dcd_s])
        if rule[0] in ["score", "casp12", "msm"]:
            cmd.append("--score")
            cmd.extend([fn.short() for fn in input_s[4]])
        if rule[0] in ["casp12"]:
            cmd.append("--qual")
            cmd.extend([fn.short() for fn in input_s[5]])
        #
        system(cmd, verbose=job.verbose)


def submit(job):
    task_s = job.get_task(METHOD, status="SUBMIT")
    if len(task_s) == 0:
        return
    #
    for index, task in task_s:
        input_s = task["input"]
        output_prefix = input_s[0]
        rule = input_s[1]
        input_json = input_s[2]
        input_dcd_s = input_s[3]
        output_pdb = task["output"][0]
        if output_pdb.status():
            continue
        #
        with input_json.open() as fp:
            options = json.load(fp)
        options["ssbond"] = []
        for line in job.ssbond:
            chain_1 = line[15]
            chain_2 = line[29]
            if chain_1 == " " and chain_2 == " ":
                line = "%sA%sA%s" % (line[:15], line[16:29], line[30:])
            options["ssbond"].append(line)
        options["rule"] = rule
        #
        job.average_d_home.chdir()
        #
        input_json = job.average_d_home.fn("%s.json" % output_prefix)
        with input_json.open("wt") as fout:
            fout.write(json.dumps(options, indent=2, default=JSONserialize))
        #
        cmd_s = []
        cmd_s.append("cd %s\n" % job.average_d_home)
        cmd = [EXEC, output_prefix, job.top_fn.short()]
        cmd.extend(["--input", input_json.short()])
        if rule[0] == "msm":
            cmd.extend(["--msm", input_s[5].short()])
        cmd.append("\\\n    ")
        cmd.append("--dcd")
        cmd.extend([fn.short() for fn in input_dcd_s])
        if rule[0] in ["score", "casp12", "msm"]:
            cmd.append("\\\n    ")
            cmd.append("--score")
            cmd.extend([fn.short() for fn in input_s[4]])
        if rule[0] in ["casp12"]:
            cmd.append("\\\n    ")
            cmd.append("--qual")
            cmd.extend([fn.short() for fn in input_s[5]])
        cmd_s.append(" ".join(cmd) + "\n")
        #
        job.write_submit_script(METHOD, index, cmd_s)


def status(job):
    pass


def main():
    arg = argparse.ArgumentParser(prog="average")
    arg.add_argument(dest="command", choices=["prep", "run"], help="exec type")
    arg.add_argument(dest="work_dir", help="work_dir, which has a JSON file")
    arg.add_argument(
        "-o", "--output", dest="output_prefix", help='output prefix, mandatory for "prep"'
    )
    arg.add_argument(
        "-i",
        "--input",
        dest="input_prod",
        type=int,
        nargs="*",
        help='input PRODs, mandatory for "prep"',
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
        if arg.input_prod is None:
            sys.exit("Error: input_prod required\n")
        #
        prep(job, arg.output_prefix, arg.input_prod)

    elif arg.command == "run":
        run(job)


if __name__ == "__main__":
    main()
