#!/usr/bin/env python

import sys
import subprocess as sp

WORK_HOME = "/home/huhlim/work/prot.str/contact/pred/trRosetta_tbm"

VERBOSE = True
PARAM_N_PROC = 16

def system(cmd, verbose=True, stdout=False, stdin=None, outfile=None, errfile=None, redirect=False):
    if type(cmd) == type(""):
        cmd = cmd.strip().split()
    if verbose and VERBOSE:
        sys.stdout.write("CMD: " + " ".join(cmd) + '\n')
    #
    if (not stdout) and (outfile is None) and (errfile is None):
        sp.call(cmd, stdin=stdin)
    else:
        if redirect:
            errfile = sp.STDOUT
        try:
            out = sp.check_output(cmd, stdin=stdin, stderr=errfile)
            if sys.version_info.major == 3:
                out = out.decode("utf8")
            if outfile is not None:
                outfile.write(out)
        except sp.CalledProcessError:
            return ''
        return out

#TBM_EXCLUDE = '%s/exclude.casp13'%WORK_HOME
TBM_EXCLUDE = None

EXEC_PSIPRED = '/home/huhlim/apps/psipred/current/run_psipred'

PARAM_SEGMENT_SIZE = 4
PARAM_PROB_CUTOFF = 0.9
PARAM_DOMAIN_CONTACT_MIN = 1000
PARAM_DOMAIN_BOUNDARY = 7
PARAM_DOMAIN_MIN_SEG = 10
PARAM_DOMAIN_MIN_RES = 30
PARAM_TBM_DOMAIN = 90.0
