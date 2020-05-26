#!/usr/bin/env python

import os
import sys
import time
import subprocess as sp

from libcasp14 import *

UPDATE_INTERVAL = 300
PARAM_BIG_TARGET = 700

def initialize_TS_server(target, fa_fn, use_hybrid=False, use_extensive=False):
    if len(target['sequence']) > PARAM_BIG_TARGET:
        send_warning(target, fa_fn)
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

def initialize_server_target(target):
    fa_fn = "%s/fa/%s.fa"%(WORK_HOME, target['target_id'])
    if not os.path.exists(fa_fn):
        with open(fa_fn, 'wt') as fout:
            fout.write(">%s\n"%(target['target_id']))
            fout.write("%s\n"%(target['sequence']))
    #
    if target['target_id'].startswith("T"):
        initialize_TS_server(target, fa_fn)
    elif target['target_id'].startswith("R"):
        initialize_TR_server(target)

def initialize_human_target(target, meta, pdb_fn, use_hybrid=False, use_extensive=False):
    EXEC = '%s/bin/casp14_refine.py'%PREFMD_HOME
    #
    cmd = [EXEC]
    cmd.append(target['target_id'])
    cmd.extend(['--input', pdb_fn])
    cmd.extend(['--dir', '%s/%s'%(WORK_HOME, meta)])
    if use_hybrid or OPTION_s[meta].get("use_hybrid", False):
        cmd.append("--hybrid")
    if use_extensive or OPTION_s[meta].get("use_extensive", False):
        cmd.append("--extensive")
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
                if META_s[meta] == 'RaptorX_TS1':
                    pdb_fn = '%s/%s.raptorX/model_1.pdb'%(TARBALL_HOME, target['target_id'])
                    if not os.path.exists(pdb_fn):
                        pdb_fn = '%s/%s/%s'%(TARBALL_HOME, target['target_id'], META_s[meta])
                else:
                    pdb_fn = '%s/%s/%s'%(TARBALL_HOME, target['target_id'], META_s[meta])
                #pdb_fn = '%s/%s/%s'%(TARBALL_HOME, target['target_id'], META_s[meta])
                if os.path.exists(pdb_fn):
                    initialize_human_target(target, meta, pdb_fn)
                    sys.stdout.write("Running %s %s\n"%(meta, target['target_id']))
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
