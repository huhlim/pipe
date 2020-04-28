#!/usr/bin/env python

import os
import sys
import path
import tempfile
import importlib
import numpy as np
from string import digits

sys.path.insert(0, '../')
from seqName import to_three_letter

class Sequence(object):
    def __init__(self, title, chain_id='', res_no=None):
        self.title = title
        self.chain_id = chain_id
        self.fasta = ''
        self.res_no = res_no
    @classmethod
    def read(cls, fa_fn):
        seq = [] ; title = ''
        with open(fa_fn) as fp:
            for line in fp:
                if line.startswith(">"):
                    title = line.strip()[1:]
                else:
                    seq.append(line.strip())
        s = cls(title)
        s.fasta = ''.join(seq)
        return s
    def write(self):
        wrt = []
        if self.res_no != None:
            if self.chain_id != '':
                wrt.append(">%s_%s_%s_%s\n"%(self.title, self.chain_id, self.res_no[0], self.res_no[1]))
            else:
                wrt.append(">%s_%s_%s\n"%(self.title, self.res_no[0], self.res_no[1]))
        else:
            if self.chain_id.strip() != '':
                wrt.append(">%s_%s\n"%(self.title, self.chain_id))
            else:
                wrt.append(">%s\n"%self.title)
        wrt.append("%s\n"%self.fasta)
        return ''.join(wrt)
    @staticmethod
    def write_SEQRES(chain_id, sequence):
        n_in_line = 13
        #
        wrt = []
        for k,resName in enumerate(sequence):
            if k%n_in_line == 0:
                wrt.append("SEQRES%4d %s %4d "%(int(k/n_in_line)+1, chain_id, len(sequence)))
            wrt.append(" %3s"%to_three_letter(resName))
            if (k+1)%n_in_line == 0:
                wrt.append("\n")
        if len(wrt) != 0 and wrt[-1] != '\n':
            wrt.append("\n")
        return wrt
    def trim_exp_tag(self):
        seq0 = self.fasta
        #
        if self.res_no is None:
            res_no = [1, len(seq0)]
        else:
            res_no = [int(self.res_no[0]), int(self.res_no[1])]
        #
        found = False ; seq = seq0
        for lib in read_EXP_TAG():
            l_seq = len(lib)
            if seq0[:l_seq] == lib:    # N-ter
                found = True
                seq = seq0[l_seq:]
                res_no[0] += l_seq
            elif seq0[-l_seq:] == lib: # C-ter
                found = True
                seq = seq0[:-l_seq]
                res_no[1] -= l_seq
            if found: break
        #
        res_no = ('%04d'%res_no[0], '%04d'%res_no[1])
        out = Sequence(self.title, chain_id=self.chain_id, res_no=res_no)
        out.fasta = seq
        return out
    def trim_signal_peptide(self, EXEC='/home/huhlim/apps/signalp/bin/signalp', options=""):
        fa_fn = tempfile.mkstemp(prefix="signalP")[1]
        with open(fa_fn, 'wt') as fout:
            fout.write(self.write())
        #
        if self.res_no is None:
            res_no = [1, len(seq0)]
        else:
            res_no = [int(self.res_no[0]), int(self.res_no[1])]
        #
        EXEC = path.Path(EXEC)
        os.environ['PATH'] = os.environ['PATH'] + ':%s'%EXEC.dirname()
        #
        cmd = []
        cmd.append('signalp')
        cmd.extend(['-fasta', fa_fn])
        cmd.extend(['-stdout', '-verbose=false'])
        cmd.extend(options.split())
        #
        sp = importlib.import_module("subprocess")
        output = sp.check_output(cmd)
        if sys.version_info.major == 3:
            output = output.decode("utf8")
        os.remove(fa_fn)
        #
        signalP = None
        for line in output.split("\n"):
            if 'CS pos: ' in line:
                signalP = int(line.split("CS pos: ")[1].split()[0].split("-")[0])
        #
        if signalP is None:
            return self
        else:
            res_no = ('%04d'%(res_no[0]+signalP), '%04d'%res_no[1])
            out = Sequence(self.title, chain_id=self.chain_id, res_no=res_no)
            out.fasta = self.fasta[signalP:]
            return out

    def __len__(self):
        return len(self.fasta)
    def __repr__(self):
        return self.fasta
    def __str__(self):
        return self.__repr__()
    def __eq__(self, othr):
        return self.fasta == othr.fasta
    def __getitem__(self, key):
        return self.fasta[key]

def read_EXP_TAG():
    lib_s = []
    #
    libdir = os.path.dirname(os.path.abspath(__file__))
    with open("%s/EXP_TAG_s"%libdir) as fp:
        for line in fp:
            lib_s.append(line.strip())
    return lib_s

read_EXP_TAG()
