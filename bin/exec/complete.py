#!/usr/bin/env python

import os
import sys
import argparse
import subprocess as sp
from string import ascii_uppercase

from genPSF import (
    AAs,
    WATERs,
    atmName_update,
    read_pdb,
    parse_patch,
    read_toppar,
    split_seg,
    write_top_cmd,
    write_pdb_cmd,
)


def write_complete_cmd():
    cmd = """bomlev -2
ic param

set tmpNIC ?NIC
coor copy comp
ic build comp
coor copy select .not. hydrogen end
hbuild atom cdie eps 80.0 cutnb 10.0 ctofnb 7.5 ctonnb 6.5 shift vshift bygr
open unit 10 write form name "out.pdb"
write coor pdb unit 10
*
close unit 10"""

    return cmd


def update_pdb(tmp_pdb, out_pdb, seg_s, ssbond_s):
    chain_id = {}
    for segName, (segType, lines, _) in seg_s.items():
        chain_id[segName] = {}
        for line in lines:
            resNo = line[22:27].strip()
            chain_id[segName][resNo] = line[21]
    #
    wrt = []
    #
    SSBOND = "SSBOND  %2d CYS %s %5s   CYS %s %5s\n"
    for disu_no, ssbond in enumerate(ssbond_s):
        chain_0 = chain_id[ssbond[0]][ssbond[1]]
        chain_1 = chain_id[ssbond[2]][ssbond[3]]
        #
        if ssbond[1][-1] in ascii_uppercase:
            cys_0_resSeq = f"{ssbond[1]:>5s}"
        else:
            cys_0_resSeq = f"{ssbond[1]:>4s} "
        if ssbond[3][-1] in ascii_uppercase:
            cys_1_resSeq = f"{ssbond[3]:>5s}"
        else:
            cys_1_resSeq = f"{ssbond[3]:>4s} "
        wrt.append(SSBOND % (disu_no + 1, chain_0, cys_0_resSeq, chain_1, cys_1_resSeq))
    #
    with open(tmp_pdb) as fp:
        for line in fp:
            if line.startswith("TER"):
                wrt.append("TER\n")
            elif line.startswith("END"):
                wrt.append("END\n")
            elif line.startswith("ATOM") or line.startswith("HETA"):
                resNo = line[22:27].strip()
                segName = line[72:76].strip()
                chain = chain_id[segName][resNo]
                wrt.append(f"{line[:21]}{chain_id[segName][resNo]}{line[22:]}")
    #
    with open(out_pdb, "wt") as fout:
        fout.writelines(wrt)


def main():
    arg = argparse.ArgumentParser(prog="complete")
    arg.add_argument(dest="init_pdb")
    arg.add_argument("-ff", "--toppar", dest="toppar", nargs="*")
    arg.add_argument("-o", "-out", "--out", dest="pdbout", default="out.pdb")
    arg.add_argument("-patch", "--patch", dest="patch_s", default=[], nargs="*")
    arg.add_argument("-blocked", "--blocked", dest="blocked", default=False, action="store_true")
    arg.add_argument("-terminal", "--terminal", dest="terminal", default=["ACE", "CT3"], nargs=2)
    arg.add_argument("--debug", dest="debug", action="store_true", default=False)
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    #
    arg.toppar = [os.path.realpath(fn) for fn in arg.toppar]
    #
    patch_s = parse_patch(arg.patch_s)
    segName_s, seg_s, disu_s, n_atoms = read_pdb(arg.init_pdb)
    tmpdir = split_seg(segName_s, seg_s, debug=arg.debug)
    pwd = os.getcwd()
    #
    cmd_s = []
    cmd_s.extend(write_top_cmd(arg.toppar))
    cmd_s.extend(
        write_pdb_cmd(
            tmpdir,
            segName_s,
            seg_s,
            disu_s,
            patch_s,
            blocked=arg.blocked,
            terminal=arg.terminal,
        )
    )
    cmd_s.append(write_complete_cmd())
    cmd_s.append("stop\n")
    #
    with open("%s/complete.cmd" % tmpdir, "wt") as fout:
        fout.writelines(cmd_s)
    #
    stdin = open("%s/complete.cmd" % tmpdir)
    stdout = open("%s/complete.log" % tmpdir, "wt")
    #
    os.chdir(tmpdir)
    if n_atoms > 360000:
        sp.call(
            os.environ["CHARMMEXEC"].split() + ["-chsize", "%d" % ((1 + int(n_atoms / 1e5)) * 1e5)],
            stdin=stdin,
            stdout=stdout,
        )
    else:
        sp.call(os.environ["CHARMMEXEC"].split(), stdin=stdin, stdout=stdout)
    os.chdir(pwd)
    #
    stdin.close()
    stdout.close()
    #
    if os.path.exists("%s/out.pdb" % tmpdir):
        update_pdb(f"{tmpdir}/out.pdb", arg.pdbout, seg_s, disu_s)
        #
        if not arg.debug:
            os.system("rm -rf %s" % tmpdir)
    else:
        sys.stdout.write("ERROR: %s\n" % tmpdir)


if __name__ == "__main__":
    main()
