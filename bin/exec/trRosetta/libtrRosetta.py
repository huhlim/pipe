#!/usr/bin/env python

import sys
import signal
import subprocess as sp

WORK_HOME = "/home/huhlim/work/prot.str/contact/pred/trRosetta_tbm"

VERBOSE = True
PARAM_N_PROC = 16

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

class GracefulExit(Exception):
    pass
def listen_signal():
    def gracefulExit(*arg):
        raise GracefulExit()
    signal.signal(signal.SIGTERM, gracefulExit)

def system(cmd, verbose=True, stdout=False, stdin=None, outfile=None, errfile=None, redirect=False):
    if type(cmd) == type(""):
        cmd = cmd.strip().split()
    if verbose:
        sys.stdout.write("CMD: " + " ".join(cmd) + '\n')
    #
    if stdout or (outfile is not None):
        STDOUT = sp.PIPE
    else:
        STDOUT = None
    if errfile is None:
        STDERR = None
    elif errfile == '/dev/null':
        STDERR = sp.DEVNULL
    elif redirect:
        STDERR = sp.STDOUT
    else:
        STDERR = errfile
    #
    listen_signal()
    #
    try:
        proc = sp.Popen(cmd, stdin=stdin, stdout=STDOUT, stderr=STDERR)
        proc.wait()
    except GracefulExit:
        proc.terminate()
        sys.exit()
    except KeyboardInterrupt:
        proc.terminate()
        sys.exit()
    #
    if STDOUT is not None:
        out = proc.stdout.read().decode("utf8")
        if outfile is not None:
            outfile.write(out)
        return out
