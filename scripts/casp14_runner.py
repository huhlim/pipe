#!/usr/bin/env python

import os
import sys
import time
import glob
import subprocess as sp

from libcasp14 import *

UPDATE_INTERVAL = 600
PARAM_BIG_TARGET = 700

def initialize_TS_server(target, fa_fn, use_hybrid=False, use_extensive=False):
    if len(target['sequence']) > PARAM_BIG_TARGET:
        send_big_target_warning(target, fa_fn)
    #
    EXEC = '%s/bin/casp14_sp.py'%PREFMD_HOME
    #
    cmd = [EXEC]
    cmd.append(target['target_id'])
    cmd.extend(['--input', fa_fn])
    cmd.extend(['--dir', '%s/%s'%(WORK_HOME, SERVER_NAME)])
    if use_hybrid or OPTION_s[SERVER_NAME][0].get("use_hybrid", False):
        cmd.append("--hybrid")
    if use_extensive or OPTION_s[SERVER_NAME][0].get("use_extensive", False):
        cmd.append("--extensive")
    sp.call(cmd)

def initialize_TR_server(target, use_hybrid=False, use_extensive=False):
    url = 'http://predictioncenter.org/download_area/CASP14/extra_experiments/%s.pdb.txt'%(target['target_id'])
    pdb_fn = '%s/refine.init/%s.pdb'%(WORK_HOME, target['target_id'])
    if not os.path.exists(pdb_fn):
        sp.call(['wget', '-q', '-c', url, '-O', pdb_fn], stderr=sp.DEVNULL)
    #
    EXEC = '%s/bin/casp14_refine.py'%PREFMD_HOME
    #
    cmd = [EXEC]
    cmd.append(target['target_id'])
    cmd.extend(['--input', pdb_fn])
    cmd.extend(['--dir', '%s/%s'%(WORK_HOME, SERVER_NAME)])
    if use_hybrid or OPTION_s[SERVER_NAME][1].get("use_hybrid", False):
        cmd.append("--hybrid")
    if use_extensive or OPTION_s[SERVER_NAME][1].get("use_extensive", False):
        cmd.append("--extensive")
    sp.call(cmd)

def initialize_TR_human(target):
    pdb_fn = '%s/refine.init/%s.pdb'%(WORK_HOME, target['target_id'])
    #
    job_home = '%s/FEIG/%s'%(WORK_HOME, target['target_id'])
    if not os.path.exists(job_home):
        os.mkdir(job_home)
    #
    init_home = '%s/init'%job_home
    if not os.path.exists(init_home):
        os.mkdir(init_home)
    if not os.path.exists("%s/init.pdb"%init_home):
        out = sp.check_output(['convpdb.pl', '-out', 'generic', pdb_fn])
        out_fn = '%s/init.pdb'%init_home
        with open(out_fn, 'wt') as fout:
            fout.write(out.decode("utf8"))
    #
    server_home = '%s/server'%job_home
    if not os.path.exists(server_home):
        os.mkdir(server_home)

    out_fn = "%s/tm.dat"%(server_home)
    if os.path.exists(out_fn):
        return
    #
    parent_id = 'T%s'%(target['target_id'].split("-")[0][1:].split("v")[0])
    server_model_s = glob.glob("%s/%s/*_TS?"%(TARBALL_HOME, parent_id))
    if len(server_model_s) == 0:
        return
    server_model_s = [fn for fn in server_model_s if not fn.split("/")[-1].startswith("server")]
    #
    cmd = []
    cmd.append("casp_eval.py")
    cmd.extend(['-r', pdb_fn])
    cmd.append("-m")
    cmd.extend(server_model_s)
    cmd.append("-simple")
    cmd.extend(['-j', '32'])
    #
    output = sp.check_output(cmd).decode("utf8").split("\n")[:-2]
    #
    wrt = [] ; dat = []
    for line in output:
        if line.startswith("#"):
            wrt.append("%s\n"%line)
        else:
            x = line.strip().split()
            gdtha = float(x[1])
            tm = float(x[3])
            dat.append((gdtha, tm, '%s\n'%line))
    dat.sort(key=lambda x: x[0], reverse=True)
    #
    pdb_fn_s = []
    for x in dat:
        wrt.append(x[2])
        if x[1] > 0.5:
            pdb_fn_s.append(x[2].strip().split()[-1])
    wrt.append("#\n")
    #
    with open(out_fn, 'wt') as fout:
        fout.writelines(wrt)
    #
    model_home = '%s/models'%server_home
    if not os.path.exists(model_home):
        os.mkdir(model_home)
    for pdb_fn in pdb_fn_s[:20]:
        sp.check_output(['cp', pdb_fn, '%s/%s.pdb'%(model_home, pdb_fn.split("/")[-1])])

def initialize_server_target(target):
    fa_fn = "%s/fa/%s.fa"%(WORK_HOME, target['target_id'])
    if not os.path.exists(fa_fn):
        with open(fa_fn, 'wt') as fout:
            fout.write(">%s\n"%(target['target_id']))
            fout.write("%s\n"%(target['sequence']))
    #
    if target['target_id'].startswith("T") or target['target_type'] == 'Private':
        initialize_TS_server(target, fa_fn)
    elif target['target_id'].startswith("R"):
        initialize_TR_server(target)
        initialize_TR_human(target)

def initialize_human_target(target_id, meta, pdb_fn, use_hybrid=False, use_extensive=False, human=None):
    EXEC = '%s/bin/casp14_refine.py'%PREFMD_HOME
    #
    cmd = [EXEC]
    cmd.append(target_id)
    cmd.extend(['--input', pdb_fn])
    cmd.extend(['--dir', '%s/%s'%(WORK_HOME, meta)])
    if use_hybrid or OPTION_s[meta].get("use_hybrid", False):
        cmd.append("--hybrid")
    if use_extensive or OPTION_s[meta].get("use_extensive", False):
        cmd.append("--extensive")
    if human is not None:
        send_action_needed_warning(target_id, meta, human)
        cmd.extend(human.split())
    sp.call(cmd)

def initialize_meta_target(target):
    EXEC = '%s/bin/casp14_meta.py'%PREFMD_HOME
    #
    cmd = [EXEC]
    cmd.append(target['target_id'])
    cmd.append("--input")
    for meta in ['FEIG-S', 'FEIG-R1', 'FEIG-R2', 'FEIG-R3']:
        job_fn = "%s/%s/%s/job.json"%(WORK_HOME, meta, target['target_id'])
        if os.path.exists(job_fn):
            cmd.append(job_fn)
    cmd.extend(['--dir', '%s/%s'%(WORK_HOME, 'FEIG')])
    sp.call(cmd)

def update_target(target, predictor):
    run_home = '%s/%s/%s'%(WORK_HOME, predictor, target['target_id'])
    final_home = '%s/final'%run_home
    #
    status = True
    for model_no in range(PARAM_N_MODEL):
        pdb_fn = '%s/model_%d.pdb'%(final_home, model_no+1)
        if not os.path.exists(pdb_fn):
            status = False
            break
    #
    job_fn = '%s/job.json'%run_home
    if os.path.exists(job_fn):
        out_s = []
        cmd = ['%s/bin/status.py'%PREFMD_HOME, 'check']
        cmd.append(job_fn)
        out_s.append(sp.check_output(cmd).decode("utf8"))
        if predictor in ['FEIG-S', 'FEIG']:
            refine_s = glob.glob("%s/refine/*/job.json"%run_home)
            refine_s.sort(key=lambda x: x.split("/")[-2])
            for refine in refine_s:
                cmd = ['%s/bin/status.py'%PREFMD_HOME, 'check', refine]
                out_s.append(sp.check_output(cmd).decode("utf8"))
        with open("%s/status.dat"%(run_home), 'wt') as fout:
            fout.write("Updated: %s\n\n"%time.ctime())
            for out in out_s:
                fout.write(out)
    #
    if not status: return False
    #
    submit_home = '%s/submit'%run_home
    if not os.path.exists(submit_home):
        os.mkdir(submit_home)
    #
    # send TS
    for model_no in range(PARAM_N_MODEL):
        pdb_fn = '%s/model_%d.pdb'%(final_home, model_no+1)
        out_fn = '%s/model_%d.pdb'%(submit_home, model_no+1)
        with open(out_fn, 'wt') as fout:
            fout.writelines(write_in_TS_format(predictor, target['target_id'], model_no, pdb_fn))
        send_prediction(predictor, target['target_id'], model_no, out_fn, mail_to=target['email'])
    #
    # send RR
    #if target['target_id'].startswith("T") and predictor == SERVER_NAME and True:
    #    npz_fn = '%s/trRosetta/%s.trRosetta.npz'%(run_home, target['target_id'])
    #    out_fn = '%s/%s.rr'%(submit_home, target['target_id'])
    #    with open(out_fn, 'wt') as fout:
    #        fout.writelines(write_in_RR_format(predictor, target['target_id'], npz_fn))
    #    send_prediction(predictor, target['target_id'], 0, out_fn, mail_to=target['email'])

    return True

def run_server(target_s):
    for target in target_s:
        if target['status'] == 'SUBMIT':
            initialize_server_target(target)
            sys.stdout.write("Running %s %s\n"%(SERVER_NAME, target['target_id']))
            #
            target['status'] = 'RUN'
            target['updated'].append("status")
        elif target['status'] == 'DONE':
            continue
        else:
            if update_target(target, SERVER_NAME):
                sys.stdout.write("Finalizing %s %s\n"%(SERVER_NAME, target['target_id']))
                target['status'] = 'DONE'
                target['updated'].append("status")

def run_human(target_s):
    for target in target_s:
        if not target['target_id'].startswith("T"):
            continue
        if target['target_type'] != 'HUMAN':
            continue
        #
        for meta in META_s:
            if META_s[meta] is None:
                continue
            #
            status_name = 'status_%s'%(meta.split("-")[1])
            if target[status_name] is None:
                multiple_domain = False ; pdb_fn_s = None
                #
                if META_s[meta] == 'RaptorX_TS1':
                    pdb_fn = '%s/%s.raptorX/model_1.pdb'%(TARBALL_HOME, target['target_id'])
                    if os.path.exists(pdb_fn):
                        pdb_fn_s = glob.glob('%s/%s.raptorX/model_1d*.pdb'%(TARBALL_HOME, target['target_id']))
                        if len(pdb_fn_s) > 0:
                            multiple_domain = True
                            pdb_fn_s.sort(key=lambda x: int(x.split("_1d")[-1][:-4]))
                    else:
                        pdb_fn = '%s/%s/%s'%(TARBALL_HOME, target['target_id'], META_s[meta])
                else:
                    pdb_fn = '%s/%s/%s'%(TARBALL_HOME, target['target_id'], META_s[meta])
                #
                if os.path.exists(pdb_fn):
                    if not multiple_domain:
                        initialize_human_target(target['target_id'], meta, pdb_fn, human=target['human'])
                        sys.stdout.write("Running %s %s\n"%(meta, target['target_id']))
                    else:
                        send_raptorX_multiple_domain_warning(target, pdb_fn_s)
                        for i,pdb_fn in enumerate(pdb_fn_s):
                            target_domain_id = '%sd%d'%(target['target_id'], i+1)
                            initialize_human_target(target_domain_id, meta, pdb_fn, human=target['human'])
                            sys.stdout.write("Running %s %s\n"%(meta, target_domain_id))
                    #
                    target[status_name] = 'RUN'
                    target['updated'].append(status_name)
            elif target[status_name] == 'DONE':
                continue
            else:
                if update_target(target, meta):
                    sys.stdout.write("Finalizing %s %s\n"%(meta, target['target_id']))
                    target[status_name] = 'DONE'
                    target['updated'].append(status_name)
        #
        if target['status_H'] is None:
            human_status = [target['status%s'%meta] == 'DONE' for meta in ['', '_R1', '_R2', '_R3']]
            if False in human_status:
                continue
            initialize_meta_target(target)
            sys.stdout.write("Running %s %s\n"%('FEIG', target['target_id']))
            #
            target['status_H'] = 'RUN'
            target['updated'].append('status_H')
        elif target['status_H'] == 'DONE':
            continue
        else:
            if update_target(target, 'FEIG'):
                sys.stdout.write("Finalizing %s %s\n"%('FEIG', target['target_id']))
                target['status_H'] = 'DONE'
                target['updated'].append('status_H')

def run():
    db, target_s = get_from_db()
    #
    get_tarball(target_s)
    get_raptorX(target_s)
    share_raptorX(target_s)
    #
    run_server(target_s)
    run_human(target_s)
    #
    update_db(db, target_s)

def main():
    while True:
        try:
            sys.stdout.write("TIME: %s\n"%time.ctime())
            run()
        except Exception as error:
            print (error)
            pass
        #
        time.sleep(UPDATE_INTERVAL)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit()
