#!/usr/bin/env python

import os
import sys
import argparse
import subprocess as sp
from tempfile import mkdtemp
from string import ascii_uppercase

import seqName
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

WORK_HOME = os.getenv("PIPE_HOME")
assert WORK_HOME is not None
sys.path.insert(0, "%s/bin" % WORK_HOME)
from libcommon import system


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


def run_scwrl(pdb_fn, tmpdir):
    input_pdb = []
    sequence = []
    with open(pdb_fn) as fp:
        for line in fp:
            if "CD  ILE" in line:
                line = line.replace("CD  ILE", "CD1 ILE")
            input_pdb.append(line)
            #
            if not line.startswith("ATOM"):
                continue
            atmName = line[12:16].strip()
            if atmName == "CA":
                resName = line[17:20]
                sequence.append(seqName.to_one_letter(resName).lower())
    sequence = "".join(sequence)
    #
    with open(f"{tmpdir}/input.pdb", "wt") as fout:
        fout.writelines(input_pdb)
    with open(f"{tmpdir}/input.fa", "wt") as fout:
        fout.write(sequence)
    #
    cmd = ["scwrl4", "-h"]
    cmd.extend(["-i", f"{tmpdir}/input.pdb"])
    cmd.extend(["-o", f"{tmpdir}/scwrl.pdb"])
    cmd.extend(["-s", f"{tmpdir}/input.fa"])
    output = system(cmd, stdout=True, verbose=False)
    #
    updated = []
    for line in output.split("\n"):
        if line.startswith("Incomplete"):
            x = line.strip().split("at")[1].split("in")[0].strip()
            updated.append((x[0], x[1:].strip()))
    #
    chain_id_prev = None
    resSeq_prev = None
    required = []
    status = []
    for line in input_pdb:
        if not line.startswith("ATOM"):
            continue
        #
        atmName = line[12:16].strip()
        resName = line[17:20]
        chain_id = line[21]
        resSeq = line[22:26].strip()
        #
        if resSeq != resSeq_prev:
            if resSeq_prev is not None:
                if False in status:
                    if (chain_id_prev, resSeq_prev) not in updated:
                        updated.append((chain_id_prev, resSeq_prev))
            #
            chain_id_prev = chain_id
            resSeq_prev = resSeq
            required = {
                "ARG": ["CZ", "NH1", "NH2"],
                "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
                "ASN": ["OD1", "ND2"],
                "GLN": ["OE1", "NE2"],
                "THR": ["OG1", "CG2"],
                "HIS": ["ND1", "CE1", "CD2", "NE2"],
            }.get(resName, [])
            status = [False for _ in required]
        if atmName in required:
            status[required.index(atmName)] = True
    #
    scwrl_pdb = {}
    with open(f"{tmpdir}/scwrl.pdb") as fp:
        for line in fp:
            if not line.startswith("ATOM"):
                continue
            #
            atmName = line[12:16].strip()
            if atmName in ["N", "CA", "C", "O"]:
                continue
            #
            chain_id = line[21]
            resSeq = line[22:26].strip()
            key = (chain_id, resSeq)
            if (chain_id, resSeq) in updated:
                if key not in scwrl_pdb:
                    scwrl_pdb[key] = []
                scwrl_pdb[key].append(line)
    #
    output = {}
    segName_s = {}
    remark = 0
    for line in input_pdb:
        if not line.startswith("ATOM"):
            remark += 1
            output[remark] = [line]
            continue
        #
        chain_id = line[21]
        resSeq = line[22:26].strip()
        key = (chain_id, resSeq)
        if key not in output:
            output[key] = []
        #
        atmName = line[12:16].strip()
        if atmName[0] == "H":
            continue
        #
        if key in updated:
            if atmName in ["N", "CA", "C", "O"]:
                segName_s[key] = line[72:76]
                output[key].append(line)
        else:
            output[key].append(line)
    for key in scwrl_pdb:
        for line in scwrl_pdb[key]:
            line = f"{line[:72]}{segName_s[key]}{line[76:]}"
            output[key].append(line)
        #
        # output[key].extend(scwrl_pdb[key])
    #
    output_pdb = []
    for key in output:
        output_pdb.extend(output[key])

    for line in output_pdb:
        if "CD1 ILE" in line:
            line = line.replace("CD1 ILE", "CD  ILE")
    #
    return output_pdb


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
    arg.add_argument("--scwrl", dest="use_scwrl", action="store_true", default=False)
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
    if arg.use_scwrl:
        if arg.debug:
            tmpdir = mkdtemp(prefix="charmm.", dir=".")
        else:
            tmpdir = mkdtemp(prefix="charmm.")
        tmpdir = os.path.abspath(tmpdir)
        #
        output = run_scwrl(arg.init_pdb, tmpdir)
        init_pdb = f"{tmpdir}/input+scwrl.pdb"
        with open(init_pdb, "wt") as fout:
            fout.writelines(output)
    else:
        init_pdb = arg.init_pdb
        tmpdir = None
    #
    segName_s, seg_s, disu_s, n_atoms = read_pdb(init_pdb)
    tmpdir = split_seg(segName_s, seg_s, tmpdir=tmpdir, debug=arg.debug)
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
