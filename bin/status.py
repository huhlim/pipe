#!/usr/bin/env python

import os
import sys
import path
import json
import time
import argparse

from libcommon import *

METHODs = ['trRosetta', 'hybrid', 'locPREFMD', 'equil', 'prod', 'score', 'average', 'qa']
INPUTs = {}
INPUTs['trRosetta'] = 2
INPUTs['hybrid'] = 2
INPUTs['locPREFMD'] = 0
INPUTs['equil'] = 0
INPUTs['prod'] = 0
INPUTs['score'] = 0
INPUTs['average'] = 0
INPUTs['qa'] = 0

def get_prod_info(run_home):
    input_json_fn = run_home.fn("input.json")
    if not input_json_fn.status():
        return 0.0, None, 0.0
    with input_json_fn.open() as fp:
        options = json.load(fp)
    dyntstep = options['md']['dyntstep']
    dynsteps = options['md']['dynsteps']
    n_iter = options['md']['iter']
    iter_simul_time = dynsteps * dyntstep / 1000.
    total_simul_time = n_iter * iter_simul_time
    #
    log_fn_s = run_home.glob("r*.*.log")
    log_fn_s.sort(key=lambda x: int(x.split('.')[-2]))
    #
    def read_log(fn):
        t_curr = 0.0 ; speed = 0.0
        with fn.open() as fp:
            for line in fp:
                if line.startswith("#") or line.startswith("ELAP"):
                    continue
                x = line.strip().split()
                try:
                    t_curr = float(x[2])
                    speed = float(x[7])
                except:
                    pass
        return t_curr / 1000., speed
    #
    if len(log_fn_s) == 0:
        return 0.0, total_simul_time, 0.0
    else:
        t_curr, speed = read_log(log_fn_s[-1])
        if speed == 0.0 and len(log_fn_s) > 1:
            _, speed = read_log(log_fn_s[-2])
        t_curr += iter_simul_time * (len(log_fn_s)-1)
        return t_curr, total_simul_time, speed

def check_status(job):
    for method in METHODs:
        new_line = False
        #
        task_s = job.get_task(method)
        if method == 'prod':
            t_now = time.time()
            for _,task in task_s:
                status = task['resource'][0]
                host = task['resource'][1]
                run_home = task['input'][0]
                #
                t_curr, t_target, speed = get_prod_info(run_home)
                if t_curr == 0.0:
                    status = 'WAIT'
                    progress = 0.0
                else:
                    if status == 'SUBMITTED':
                        status = 'RUN'
                    progress = min(1.0, t_curr / t_target) * 100.0
                if speed > 0.0:
                    t_left = (t_target-t_curr) / speed * 86400.
                    t_est = time.ctime(t_now + t_left)
                else:
                    t_left = None
                    t_est = None
                #
                wrt = []
                wrt.append("%-10s"%job.title)
                wrt.append('%-10s'%method)
                wrt.append("%-6s"%status)
                wrt.append("%-12s"%host)
                wrt.append("%s"%(run_home.short()))
                wrt.append("%5.1f [%%]"%progress)
                wrt.append("%6.1f [ns/day]"%speed)
                if t_left is not None and status == 'RUN':
                    wrt.append("%6.2f hrs_left"%(t_left/3600.))
                    wrt.append(t_est)
                sys.stdout.write("  ".join(wrt) + '\n')
                new_line = True
        else:
            for _,task in task_s:
                status = task['resource'][0]
                host = task['resource'][1]
                if method != 'average':
                    input = task['input'][INPUTs[method]].short()
                else:
                    input = task['input'][INPUTs[method]]
                #
                wrt = []
                wrt.append("%-10s"%job.title)
                wrt.append('%-10s'%method)
                wrt.append("%-6s"%status)
                wrt.append("%-12s"%host)
                wrt.append("%s"%input)
                sys.stdout.write("  ".join(wrt) + '\n')
                new_line = True
        if new_line:
            sys.stdout.write("#\n")

def check_resource(host_s, json_job_s):
    resource_s = {}
    for host in host_s:
        resource_s[host] = {}
        resource_s[host]['cpu'] = None
        if host_s[host][0]:
            resource_s[host]['gpu'] = {}
            for gpu_id in host_s[host][1]:
                resource_s[host]['gpu'][gpu_id] = None
    #
    for json_job in json_job_s:
        if not json_job.status():
            continue
        #
        job = Job.from_json(json_job)
        for method in job.task:
            task_s = job.get_task(method, status='RUN')
            for index, task in task_s:
                host = task['resource'][1]
                is_gpu = task['resource'][2]
                if method != 'average':
                    input = task['input'][INPUTs[method]].short()
                else:
                    input = task['input'][INPUTs[method]]
                if is_gpu:
                    host_name = host.split("/")[0]
                    gpu_id = int(host.split("/")[1])
                    resource_s[host_name]['gpu'][gpu_id] = (job.title, method, input)
                else:
                    resource_s[host]['cpu'] = (job.title, method, input)
    #
    for host in sorted(host_s):
        if 'gpu' not in resource_s[host]: continue
        #
        for gpu_id in sorted(resource_s[host]['gpu']):
            if resource_s[host]['gpu'][gpu_id] is None:
                wrt = ''
            else:
                wrt = '%-10s  %-10s  %s'%resource_s[host]['gpu'][gpu_id]
            sys.stdout.write("%-10s GPU %d   %s\n"%(host, gpu_id, wrt))
    #
    sys.stdout.write("#\n")
    #
    for host in sorted(host_s):
        if resource_s[host]['cpu'] is None:
            wrt = ''
        else:
            wrt = '%-10s  %-10s  %s'%resource_s[host]['cpu']
        sys.stdout.write("%-10s CPU     %s\n"%(host, wrt))



def main():
    arg = argparse.ArgumentParser(prog='PREFMD.status')
    sub = arg.add_subparsers(dest='method')
    #
    sub_list = sub.add_parser("list")
    #
    sub_check = sub.add_parser("check")
    sub_check.add_argument("json_job_s", nargs='*')

    sub_update = sub.add_parser("update")

    sub_append = sub.add_parser("append")
    sub_append.add_argument("json_job_s", nargs='+')

    sub_remove = sub.add_parser("remove")
    sub_remove.add_argument("json_job_s", nargs='+')

    sub_resource = sub.add_parser("resource")

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    if arg.method == 'list':
        if JOBs_json.status():
            with JOBs_json.open() as fp:
                json_job_s = [path.Path(fn) for fn in json.load(fp)]
        else:
            return
        for json_job in json_job_s:
            sys.stdout.write("%s\n"%json_job)

    elif arg.method == 'check':
        json_job_s = arg.json_job_s
        if len(json_job_s) == 0:
            if JOBs_json.status():
                with JOBs_json.open() as fp:
                    json_job_s = [path.Path(fn) for fn in json.load(fp)]

        for json_job_fn in json_job_s:
            if json_job_fn.endswith(".json"):
                json_job_fn = path.Path(json_job_fn)
            else:
                json_job_fn = path.Dir(json_job_fn).fn("job.json")
            #
            job = Job.from_json(json_job_fn)
            check_status(job)

    elif arg.method == 'update':
        return

    elif arg.method == 'append':
        if JOBs_json.status():
            with JOBs_json.open() as fp:
                job_s = [fn for fn in json.load(fp)]
        else:
            job_s = []

        job_s += [path.Path(fn).path() for fn in arg.json_job_s]
        #
        with JOBs_json.open('wt') as fout:
            fout.write(json.dumps(job_s, indent=2))

    elif arg.method == 'remove':
        json_job_s = [path.Path(fn) for fn in arg.json_job_s]
        if JOBs_json.status():
            with JOBs_json.open() as fp:
                job_s = []
                for fn in json.load(fp):
                    json_job = path.Path(fn)
                    if json_job not in json_job_s:
                        job_s.append(fn)
            with JOBs_json.open('wt') as fout:
                fout.write(json.dumps(job_s, indent=2))

    elif arg.method == 'resource':
        if HOSTs_json.status():
            with HOSTs_json.open() as fp:
                host_s = json.load(fp)
        else:
            return
        if JOBs_json.status():
            with JOBs_json.open() as fp:
                json_job_s = [path.Path(fn) for fn in json.load(fp)]
        else:
            json_job_s = []
        #
        check_resource(host_s, json_job_s)

if __name__ == '__main__':
    main()
