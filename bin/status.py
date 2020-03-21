#!/usr/bin/env python

import os
import sys
import path
import json
import time

from libcommon import *

METHODs = ['hybrid', 'locPREFMD', 'equil', 'prod', 'score', 'average', 'qa']
INPUTs = {}
INPUTs['hybrid'] = 2
INPUTs['locPREFMD'] = 0
INPUTs['equil'] = 0
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

def run(job):
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
                wrt.append("%-10s"%host)
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
                wrt.append("%-10s"%host)
                wrt.append("%s"%input)
                sys.stdout.write("  ".join(wrt) + '\n')
                new_line = True
        if new_line:
            sys.stdout.write("#\n")

def main():
    json_job_s = sys.argv[1:]
    for json_job_fn in json_job_s:
        if json_job_fn.endswith(".json"):
            json_job_fn = path.Path(json_job_fn)
        else:
            json_job_fn = path.Dir(json_job_fn).fn("job.json")
        #
        job = Job.from_json(json_job_fn)
        run(job)

if __name__ == '__main__':
    main()
