#!/usr/bin/env python

import os
import sys
import datetime
import pymysql
import requests
import numpy as np
import subprocess as sp

PREFMD_HOME = os.getenv("PREFMD_HOME")
WORK_HOME = '/green/s2/huhlim/work/casp14'
PARAM_N_MODEL = 5

PARAM_CASP_EMAIL = 'models@predictioncenter.org'
PARAM_MY_EMAIL = 'huhlim@gmail.com'
TARBALL_HOME = '%s/servers'%WORK_HOME

SERVER_NAME = 'FEIG-S'

CODE_s = {}
CODE_s['FEIG']    = '8100-2136-8593'
CODE_s['FEIG-S']  = '4895-5981-8500'
CODE_s['FEIG-R1'] = '2545-9941-4948'
CODE_s['FEIG-R2'] = '6341-4498-6818'
CODE_s['FEIG-R3'] = '3984-5248-0113'
CODE_s['FEIG-R4'] = '7866-2442-6941'

METHOD_s = {}
METHOD_s['FEIG']    = ['FeigLab structure prediction', 'FeigLab refinement']
METHOD_s['FEIG-S']  = ['FeigLab structure prediction', 'FeigLab refinement']
METHOD_s['FEIG-R1'] = ['RaptorX + PREFMD', '']
METHOD_s['FEIG-R2'] = ['Zhang-Server + PREFMD', '']
METHOD_s['FEIG-R3'] = ['BAKER-ROSETTASERVER + PREFMD', '']
METHOD_s['FEIG-R4'] = ['', '']

META_s = {}
META_s['FEIG-R1'] = 'RaptorX_TS1'
META_s['FEIG-R2'] = 'Zhang-Server_TS1'
META_s['FEIG-R3'] = 'BAKER-ROSETTASERVER_TS1'
META_s['FEIG-R4'] = None

OPTION_s = {}
OPTION_s['FEIG']    = {}
OPTION_s['FEIG-S']  =[{"use_hybrid": True}, \
                      {"use_hybrid": True}]
OPTION_s['FEIG-R1'] = {}
OPTION_s['FEIG-R2'] = {}
OPTION_s['FEIG-R3'] = {}
OPTION_s['FEIG-R4'] = {}

def send_prediction(predictor, target_id, model_no, fn, mail_to=PARAM_CASP_EMAIL):
    cmd = []
    cmd.append("mail")
    cmd.append("-s")
    cmd.append("%s MODEL%d - %s"%(target_id, model_no+1, predictor))
    cmd.append("-r")
    cmd.append(PARAM_MY_EMAIL)
    cmd.append(mail_to)
    #
    with open(fn) as content:
        sp.check_output(cmd, stdin=content)

def send_warning(target, fn):
    cmd = []
    cmd.append("mail")
    cmd.append("-s")
    cmd.append("A big CASP14 target has submitted: %s (%d)"%(target['target_id'], len(target['sequence'])))
    cmd.append("-r")
    cmd.append(PARAM_MY_EMAIL)
    cmd.append(PARAM_MY_EMAIL)
    #
    with open(fn) as content:
        sp.check_output(cmd, stdin=content)

def write_in_TS_format(predictor, target_id, model_no, pdb_fn):
    pdb = []
    with open(pdb_fn) as fp:
        for line in fp:
            if line.startswith("ATOM"):
                pdb.append(line)
            elif line.startswith("TER"):
                pdb.append("TER\n")
    #
    wrt = []
    wrt.append("PFRMAT TS\n")
    wrt.append("TARGET %s\n"%target_id)
    wrt.append("AUTHOR %s\n"%CODE_s[predictor])
    if target_id.startswith("R"):
        wrt.append("METHOD %s\n"%METHOD_s[predictor][1])
    else:
        wrt.append("METHOD %s\n"%METHOD_s[predictor][0])
    wrt.append("MODEL %d\n"%(model_no+1))
    wrt.append("PARENT N/A\n")
    wrt.extend(pdb)
    if wrt[-1].strip() != 'TER':
        wrt.append("TER\n")
    wrt.append("END\n")
    return wrt

def write_in_RR_format(predictor, target_id, npz_fn):
    PARAM_topL = 10 
    distIndex = np.arange(36, dtype=float) * 0.5 + 2.25
    distIndex = ((distIndex-2.0) / 2.0).astype(np.int)
    distIndex = np.concatenate([[9], distIndex])

    prob = np.load(npz_fn)['dist']
    l_seq = prob.shape[0]
    prob_10 = np.zeros((l_seq,l_seq,10), dtype=np.float)
    for i in range(l_seq-2):
        for j in range(i+2,l_seq):
            np.add.at(prob_10[i,j], distIndex, prob[i,j])

    prob_2 = np.sum(prob_10[:,:,:3], axis=-1)

    n_out = PARAM_topL * l_seq

    prob_cut = np.sort(np.concatenate(prob_2))[::-1][n_out]
    prob_cut = np.maximum(prob_cut, 0.001)
    #
    pair_s = np.where(prob_2 >= prob_cut)
    #
    wrt_s = []
    wrt_s.append("PFRMAT RR\n")
    wrt_s.append("TARGET %s\n"%target_id)
    wrt_s.append("AUTHOR %s\n"%CODE_s[predictor])
    wrt_s.append("METHOD trRosetta + FeigLab modification\n")
    wrt_s.append("RMODE 1\n")
    wrt_s.append("MODEL 1\n")
    for i,j in zip(pair_s[0], pair_s[1]):
        wrt = []
        wrt.append("%d %d"%(i+1, j+1))
        wrt.append("%5.3f"%prob_2[i,j])
        wrt.append(' '.join(['%5.3f'%p for p in prob_10[i,j]]))
        wrt = '  '.join(wrt)
        wrt_s.append(wrt + '\n')
    wrt_s.append("END\n")
    return wrt_s

def get_from_db():
    db = pymysql.connect(host='blue7', user='casp14', password='CA9ee04cf01SP14!', db='casp14',
                cursorclass=pymysql.cursors.DictCursor)
    #
    with db.cursor() as cursor:
        sql = 'SELECT * from targets'
        cursor.execute(sql)
        target_s = cursor.fetchall()
    #
    for target in target_s:
        target['updated'] = []
    #
    parse_target_list(target_s)
    #
    return db, target_s

def update_db(db, target_s):
    need_commit = False
    with db.cursor() as cursor:
        for target in target_s:
            target['updated'] = list(set(target['updated']))
            if len(target['updated']) == 0:
                continue
            for column in target['updated']:
                sql = 'UPDATE targets SET'
                sql += ' %s = "%s"'%(column, target[column])
                sql += ' where id = %d'%(target['id'])
                cursor.execute(sql)
                need_commit = True
    #
    if need_commit:
        db.commit()
    db.close()

def parse_target_list(target_s):
    know_target_types = True
    for target in target_s:
        if target['target_type'] is None:
            know_target_types = False
            break
    if know_target_types:
        return
    #
    target_type_s = {}
    csv = requests.get("http://predictioncenter.org/casp14/targetlist.cgi?type=csv").text
    for line in csv.split("\n"):
        if line.startswith("Target;"):
            continue
        x = line.strip().split(";")
        if len(x) < 2: continue
        #
        target_id = x[0]
        if x[1].startswith("All"):
            target_type_s[target_id] = 'HUMAN'
        else:
            target_type_s[target_id] = 'SERVER'
    #
    for target in target_s:
        if target['target_type'] is not None:
            continue
        if target['target_id'] in target_type_s:
            target['target_type'] = target_type_s[target['target_id']]
        else:
            target['target_type'] = "CASP13"
        target['updated'].append("target_type")

def get_tarball(target_s):
    cwd = os.getcwd()
    now = datetime.datetime.now()
    #
    for target in target_s:
        if target['target_type'] not in ['HUMAN', 'SERVER']:
            continue
        #
        run_home = '%s/%s'%(TARBALL_HOME, target['target_id'])
        if not os.path.exists(run_home):
            os.mkdir(run_home)
        #
        submit_time = target['submit_time']
        #
        for tarball_days,tarball_type in zip([7, 9], ['stage2.3D', '3D']):
            tarball_time = submit_time + datetime.timedelta(days=tarball_days)
            tarball_time = tarball_time.replace(hour=15, minute=0, second=0, microsecond=0)
            if now < tarball_time:
                continue
            #
            tgz_fn = '%s.%s.srv.tar.gz'%(target['target_id'], tarball_type)
            if os.path.exists('%s/%s'%(run_home, tgz_fn)):
                continue
            #
            url = 'http://predictioncenter.org/download_area/CASP14/server_predictions/%s'%tgz_fn
            sp.call(['wget', '-q', '-c', url, '-O', '%s/%s'%(run_home, tgz_fn)], stderr=sp.DEVNULL)
            if not os.path.exists('%s/%s'%(run_home, tgz_fn)):
                continue
            if os.path.getsize('%s/%s'%(run_home, tgz_fn)) < 1024:
                os.remove('%s/%s'%(run_home, tgz_fn))
                continue
            #
            os.chdir(run_home)
            sp.call(['tar', 'xzf', tgz_fn, '-C', '..'])
    #
    os.chdir(cwd)
