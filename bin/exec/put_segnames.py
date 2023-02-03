#!/usr/bin/env python

import os
import sys

from libmd import put_segnames


def main():
    if len(sys.argv) < 3:
        sys.exit(f"usage: {__file__} [in_pdb] [out_pdb] (chain_break)")
    #
    in_pdb = sys.argv[1]
    out_pdb = sys.argv[2]
    if len(sys.argv) > 3:
        chain_break_distance = float(sys.argv[3])
    else:
        chain_break_distance = 2.0
    #
    out, segName_s = put_segnames(in_pdb, CHAIN_BREAKs=chain_break_distance)
    with open(out_pdb, "wt") as fout:
        fout.writelines(out)


if __name__ == "__main__":
    main()
