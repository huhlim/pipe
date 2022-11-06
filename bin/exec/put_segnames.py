#!/usr/bin/env python

import os
import sys

from libmd import put_segnames


def main():
    if len(sys.argv) < 3:
        sys.exit(f"usage: {__file__} [in_pdb] [out_pdb]")
    #
    in_pdb = sys.argv[1]
    out_pdb = sys.argv[2]
    #
    out = put_segnames(in_pdb)
    with open(out_pdb, "wt") as fout:
        fout.writelines(out)


if __name__ == "__main__":
    main()
