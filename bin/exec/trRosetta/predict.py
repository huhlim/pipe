#!/usr/bin/env python

import os
import sys
import copy
import path
import mdtraj
import numpy as np
import networkx as nx

from libseq import Sequence

from libtrRosetta import *

from run_tbm import run as run_tbm
from run_tbm import read_tbm
from build_msa import run as build_msa

from run_trRosetta import read_trRosetta, get_contact_trRosetta, FEATUREs
from run_trRosetta import merge as merge_trRosetta
from run_trRosetta import run as run_trRosetta

from get_domain import run as get_domain
from tbm_to_contact import run as tbm_to_contact

from build_model import run as build_model 

class Job(object):
    def __init__(self, run_home, fa_fn):
        self.fa_fn0 = path.Path(fa_fn)
        self.title = path.name(fa_fn)
        self.run_home = path.Dir(run_home, build=True)
        self.domain_s = BreadthFirstSearch()
    def __repr__(self):
        return self.title
    def initialize_run(self, forced=False):
        self.DONE = self.fn("DONE")
        self.LOCK = self.fn("LOCK")
        #
        if not forced:
            if self.DONE.status(): return False
            if self.LOCK.status(): return False
        #
        with self.LOCK.open('wt') as fout:
            fout.write("#")
        self.DOMAIN = self.fn("domain_s")
        return True
    def finalize_run(self):
        with self.DONE.open('wt') as fout:
            fout.write("#")
        self.LOCK.remove()
    def fn(self, X=None, suffix=None):
        if X is not None:
            return self.run_home.fn(X)
        elif suffix is not None:
            return self.run_home.fn("%s.%s"%(self.title, suffix))
    def trim_artificial_sequence(self, run=True):
        seq0 = Sequence.read(self.fa_fn0.short())
        self.seq0 = seq0
        #
        if run:
            self.fa_fn = self.fn(suffix='fa')
            if not self.fa_fn.status():
                seq1 = seq0.trim_exp_tag()
                seq2 = seq1.trim_signal_peptide()
                #
                with self.fa_fn.open("wt") as fout:
                    fout.write(seq2.write())
            else:
                seq2 = Sequence.read(self.fa_fn.short())
        else:
            self.fa_fn = self.fa_fn0
            seq2 = seq0
        self.seq  = seq2
        #
        self.searched = np.zeros(len(self.seq0), dtype=bool)
        n_ter = self.seq0.fasta.find(self.seq.fasta)
        self.searched[n_ter:n_ter+len(self.seq)] = True
    def run_psipred(self):
        def read_psipred(psipred_fn):
            ss = []
            with psipred_fn.open() as fp:
                for line in fp:
                    if line.startswith("#"):
                        continue
                    line = line.strip()
                    if len(line) == 0:
                        continue
                    ss.append(line.split()[2])
            ss = ''.join(ss)
            return ss
        #
        psipred_home = self.run_home.subdir("ss", build=True)
        psipred_fn = psipred_home.fn("%s.ss2"%self.title)
        #
        if not psipred_fn.status():
            psipred_home.chdir()
            #
            cmd = []
            cmd.append(EXEC_PSIPRED)
            cmd.append(self.fa_fn.short())
            system(cmd)
        #
        self.run_home.chdir()
        self.ss = read_psipred(psipred_fn)

    def add_new_domain(self, mask, parent=None):
        is_new = True
        for domain in self.domain_s:
            if np.all(mask == domain.mask):
                is_new = False ; break
        if not is_new: return
        #
        nth_domain_search = len(self.domain_s)
        if parent is None:
            name = 'domain_%d_%d'%(0, nth_domain_search)
        else:
            name = 'domain_%d_%d'%(parent.search_no, nth_domain_search)
        #
        domain = Domain(nth_domain_search, self.title, name, self.seq, mask, parent=parent)
        domain.domain_home(run_home=self.run_home)
        self.domain_s.append(domain)
        self.write_domain_info()
    def write_domain_info(self):
        with self.DOMAIN.open('wt') as fout:
            for domain in self.domain_s:
                fout.write('%s\n'%domain)
    def read_domain_info(self):
        if not self.DOMAIN.status():
            return
        with self.DOMAIN.open() as fp:
            for line in fp:
                name, mask_str = line.strip().split()
                mask = np.array([X == 'X' for X in mask_str], dtype=bool)
                #
                x = name.split("_")
                parent_id = int(x[1])
                nth_domain_search = int(x[2])
                if nth_domain_search == 0:
                    parent = None
                else:
                    parent = self.domain_s[parent_id]
                    if parent.fn("trRosetta/%s.npz"%self.title).status():
                        parent.searched = True
                #
                domain = Domain(nth_domain_search, self.title, name, self.seq, mask, parent=parent)
                domain.domain_home(run_home=self.run_home)
                self.domain_s.append(domain)
            #
            for domain in self.domain_s:
                if domain.searched: continue
                if domain.fn("trRosetta/%s.npz"%self.title).status():
                    domain.searched = True

class BreadthFirstSearch(object):
    def __init__(self):
        self.domain_s = {}
        self.graph = nx.DiGraph()
    def __len__(self):
        return len(self.domain_s)
    def __iter__(self):
        self.key_s = sorted(self.domain_s.keys())
        self.iter_index = -1
        self.n_iter = len(self.key_s)
        return self
    def __next__(self):
        self.iter_index += 1
        if self.iter_index == self.n_iter:
            raise StopIteration
        key = self.key_s[self.iter_index]
        return self.domain_s[key]
    def __getitem__(self, index):
        return self.domain_s[index]
    def append(self, domain):
        if domain.parent is None:
            self.graph.add_node(domain.search_no)
        else:
            self.graph.add_edge(domain.parent.search_no, domain.search_no)
        self.domain_s[domain.search_no] = domain
    def get_unsearched(self):
        if not self.domain_s[0].searched:
            return self.domain_s[0]
        for search_no, _ in nx.bfs_predecessors(self.graph, 0):
            if not self.domain_s[search_no].searched:
                return self.domain_s[search_no]
        return None
    def has_unsearched(self):
        if not self.domain_s[0].searched:
            return True
        for search_no, _ in nx.bfs_predecessors(self.graph, 0):
            if not self.domain_s[search_no].searched:
                return True
        return False
    def get_splitted(self):
        if self.domain_s[0].is_tbm:
            return [0]
        #
        def has_tbm_domain(self, parent):
            for search_no, _ in nx.bfs_predecessors(self.graph, parent):
                if self.domain_s[search_no].is_tbm:
                    return True
            return False
        #
        if not has_tbm_domain(self, 0):
            return [0]
        #
        splitted = [] ; successors = []
        for search_no, _ in nx.bfs_predecessors(self.graph, 0):
            if search_no in successors: continue
            #
            has_successor = (len(next(nx.bfs_successors(self.graph, search_no))[1]) > 0)
            if self.domain_s[search_no].is_tbm:
                splitted.append(search_no)
                for x,_ in nx.bfs_predecessors(self.graph, search_no):
                    successors.append(x)
            elif not has_successor:
                splitted.append(search_no)
                for x,_ in nx.bfs_predecessors(self.graph, search_no):
                    successors.append(x)
            elif not has_tbm_domain(self, search_no):
                splitted.append(search_no)
                for x,_ in nx.bfs_predecessors(self.graph, search_no):
                    successors.append(x)
        return splitted

class Domain(object):
    def __init__(self, search_no, title, name, seq0, mask, parent=None):
        self.search_no = search_no
        self.title = title
        self.name = name
        self.mask = mask
        self.searched = False
        self.tbm = None
        self.parent = parent
        #
        self.seq = Sequence(title)
        self.seq.fasta = ''.join([x for x,m in zip(seq0, mask) if m])
    def __repr__(self):
        wrt = ['%-12s'%self.name]
        mask_str = []
        for x in self.mask:
            if x: mask_str.append("X")
            else: mask_str.append("-")
        wrt.append("".join(mask_str))
        return '  '.join(wrt)
    def domain_home(self, run_home=None):
        if hasattr(self, '_domain_home'):
            return self._domain_home
        else:
            assert run_home != None
            self._domain_home = run_home.subdir(self.name, build=True)
            return self._domain_home
    def write_fasta(self, fa_fn):
        with fa_fn.open("wt") as fout:
            fout.write(self.seq.write())
    @property
    def is_tbm(self):
        if self.tbm is not None:
            return self.tbm
        tbm_s = read_tbm(self.title, self.domain_home())
        for method in tbm_s:
            for templ in tbm_s[method]:
                if templ[2] > PARAM_TBM_DOMAIN:
                    self.tbm = True
                    return self.tbm
        self.tbm = False
        return self.tbm

    def __hash__(self):
        return self.name.__hash__()
    def __eq__(self, othr):
        if isinstance(othr, Domain):
            return self.name == othr.name
        else:
            return self.name == othr
    def fn(self, X=None, suffix=None):
        if X is not None:
            return self._domain_home.fn(X)
        elif suffix is not None:
            return self._domain_home.fn("%s.%s"%(self.title, suffix))

def hybrid_contact(npz_s):
    hybrid = {}
    for feature in FEATUREs:
        hybrid[feature] = np.zeros_like(npz_s[0][1][feature])
        l_seq = hybrid[feature][0].shape[0]
        weight_sum = np.zeros((l_seq, l_seq), dtype=np.float32)
        #
        for w,npz in npz_s:
            if isinstance(w, float):
                hybrid[feature] += w*npz[feature]
            else:
                hybrid[feature] += w[:,:,None]*npz[feature]
            weight_sum += w
        hybrid[feature] /= weight_sum[:,:,None]
    return hybrid

def run(job):
    job.trRosetta_fn = job.fn(suffix='trRosetta.npz')
    if job.trRosetta_fn.status():
        return
    if len(job.domain_s) == 0:
        job.add_new_domain(np.ones(len(job.seq), dtype=bool))
    #
    while job.domain_s.has_unsearched():
        domain = job.domain_s.get_unsearched()
        mask0 = domain.mask
        #
        fa_fn = domain.fn(suffix='fa')
        if not fa_fn.status():
            domain.write_fasta(fa_fn)
        #
        # TBM-search
        tbm_s = run_tbm(fa_fn, n_proc=PARAM_N_PROC)
        #
        # MSA generation
        msa_s = build_msa(fa_fn, n_proc=PARAM_N_PROC)
        #
        # trRosetta
        trRosetta = run_trRosetta(msa_s['trRosetta'])
        contact_trRosetta = get_contact_trRosetta(trRosetta)
        #
        # domain_detection
        domain_ss = ''.join([ss for ss,m in zip(job.ss, mask0) if m])
        domain_new_s = get_domain(domain_ss, contact_trRosetta)
        for new in domain_new_s:
            mask = copy.deepcopy(mask0)
            mask[mask] = new
            job.add_new_domain(mask, parent=domain)
        #
        domain.searched = True

    if not job.trRosetta_fn.status():
        merged = None
        for domain in job.domain_s:
            npz_fn = domain.fn("trRosetta/%s.npz"%job.title)
            if not npz_fn.status():
                sys.exit("Failed to find %s\n"%npz_fn)
            #
            npz_TBM = tbm_to_contact(job.title, domain)
            npz_trRosetta = read_trRosetta(npz_fn)
            npz = [(0.5, npz_trRosetta)] + npz_TBM
            npz = hybrid_contact(npz)

            merged = merge_trRosetta(merged, npz, domain)
        #
        out_s = {} ; l_seq = len(job.seq0)
        for feature in FEATUREs:
            out = np.zeros((l_seq, l_seq, merged[feature].shape[-1]), dtype=np.float32)
            mask = np.ix_(job.searched, job.searched)
            out[mask] = merged[feature]
            out_s[feature] = out

        np.savez(job.trRosetta_fn.short(), **out_s)
    #
    job.run_home.chdir()

def split_model(job, pdb_fn):
    model_home = job.run_home.subdir("model", build=True)
    model_home.chdir()
    #
    pdb0 = mdtraj.load(pdb_fn.short())
    resNo = np.array([atom.residue.index for atom in pdb0.top.atoms])
    #
    l_seq = len(job.seq0)
    out_fn_s = []
    for split in job.domain_s.get_splitted():
        out_fn = model_home.fn("%s_d%d.pdb"%(job.title, split))
        if out_fn.status():
            continue
        mask = np.zeros(l_seq, dtype=bool)
        mask[job.searched] = job.domain_s[split].mask
        selected = np.where(mask[resNo])[0]
        #
        pdb = pdb0.atom_slice(selected)
        pdb.save(out_fn.short())
    #
    with job.fn("model_s").open("wt") as fout:
        for out in out_fn_s:
            fout.write("%s\n"%out)

def main():
    if len(sys.argv) == 1:
        sys.stderr.write("USAGE: %s [FA]\n"%__file__)
        return
    #
    run_home = os.getcwd()
    in_fa = sys.argv[1]
    #
    job = Job(run_home, in_fa)
    if not job.initialize_run():
        return
    #
    job.trim_artificial_sequence()
    job.run_psipred()
    job.read_domain_info()
    #
    run(job)
    pdb_fn = build_model(job)
    split_model(job, pdb_fn)
    #
    job.finalize_run()
    #try:
    #except:
    #    job.LOCK.remove()

if __name__ == '__main__':
    main()
