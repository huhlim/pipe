#!/usr/bin/env python

import os
import sys
import argparse
import subprocess as sp

from libcasp14 import *

STATUS_NAMEs = {}
STATUS_NAMEs['FEIG-S'] = 'status'
STATUS_NAMEs['FEIG-R1'] = 'status_R1'
STATUS_NAMEs['FEIG-R2'] = 'status_R2'
STATUS_NAMEs['FEIG-R3'] = 'status_R3'
STATUS_NAMEs['FEIG-R4'] = 'status_R4'
STATUS_NAMEs['FEIG'] = 'status_H'

def send_FINAL_model(predictor, target_id):
    run_home = '%s/%s/%s'%(WORK_HOME, predictor, target_id)
    #
    final_home = '%s/final'%run_home
    for model_no in range(PARAM_N_MODEL):
        out_fn = '%s/model_%d.pdb'%(final_home, model_no+1)
        if not os.path.exists(out_fn):
            sys.stderr.write("Error: missing %s\n"%out_fn)
            return
    #
    submit_home = '%s/submit'%run_home
    if not os.path.exists(submit_home):
        os.mkdir(submit_home)
    #
    for model_no in range(PARAM_N_MODEL):
        out_fn = '%s/model_%d.pdb'%(final_home, model_no+1)
        pdb_fn = '%s/model_%d.pdb'%(submit_home, model_no+1)
        #
        with open(pdb_fn, 'wt') as fout:
            fout.writelines(write_in_TS_format(predictor, target_id, model_no, out_fn))
        #
        send_prediction(predictor, target_id, model_no, pdb_fn)
        sys.stdout.write("SENDING ... %s %s model_%d\n"%(predictor, target_id, model_no+1))

def main():
    arg = argparse.ArgumentParser(prog='casp14send')
    arg.add_argument(dest='predictor', help='predictor', choices=CODE_s.keys())
    arg.add_argument(dest='target_id', help='Job title')

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    #
    send_FINAL_model(arg.predictor, arg.target_id)
    #
    db, target_s = get_from_db()
    for target in target_s:
        if target['target_id'] != arg.target_id:
            continue
        status_name = STATUS_NAMEs[arg.predictor]
        target[status_name] = 'DONE'
        target['updated'].append(status_name)
    update_db(db, target_s)

if __name__ == '__main__':
    main()
