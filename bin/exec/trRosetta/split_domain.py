#!/usr/bin/env python

import os
import sys
import networkx as nx
from predict import Job

def main():
    if len(sys.argv) == 1:
        sys.stderr.write("USAGE: %s [TEST] [FA]\n"%__file__)
        return
    #
    run_home = sys.argv[1]
    in_fa = sys.argv[2]
    #
    job = Job(run_home, in_fa)
    if not job.initialize_run(forced=True):
        return
    #
    job.trim_artificial_sequence()
    job.run_psipred()
    job.read_domain_info()
    #
    for split in job.domain_s.get_splitted():
        print (job.domain_s[split].mask)

if __name__ == '__main__':
    main()

