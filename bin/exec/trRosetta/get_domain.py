#!/usr/bin/env python

import os
import sys
import copy
import path
import numpy as np
import networkx as nx
from scipy.sparse import csgraph

from libseq import Sequence

from libtrRosetta import *

class SegmentMap(object):
    def __init__(self, ss):
        self.ss = ss
        self.l_seq = len(ss)
        self.define_segments()
    def __repr__(self):
        return self.ss
    def define_segments(self):
        segment = np.zeros(len(self.ss), dtype=int)
        #
        seg_no = 0 ; l_seg = 0 ; ss_prev = None
        for i,s in enumerate(self.ss):
            if s != ss_prev:
                seg_no += 1
            elif l_seg == PARAM_SEGMENT_SIZE:
                seg_no += 1
                l_seg = 0
            l_seg += 1
            ss_prev = s
            segment[i] = seg_no
        segment -= 1
        #
        self.segment_index = segment
        self.n_segment = seg_no
        self.map0 = np.eye(self.n_segment, dtype=float)
        #
        for i in range(self.n_segment-1):
            x = np.where(self.segment_index==i)[0][-1]
            y = np.where(self.segment_index==i+1)[0][0]
            if self.ss[x] == self.ss[y] and (y-x == 1):   # in a same SSE
                self.map0[i,i+1] = 1.
                self.map0[i+1,i] = 1.

    def apply_map(self, prob, verbose=False):
        def get_domain(self, Xs, segment_i, verbose=False):
            domain_s = []
            for xi in Xs:
                out = np.zeros(len(self.ss), dtype=bool)
                for l in xi:
                    out[self.segment_index == segment_i[l]] = True
                out_ext = np.zeros_like(out)
                for i,selected in enumerate(out):
                    if not selected:
                        continue
                    r_low = max(0, i-PARAM_DOMAIN_BOUNDARY)
                    r_max = min(self.l_seq, i+PARAM_DOMAIN_BOUNDARY)
                    out_ext[r_low:r_max] = True
                domain_s.append(out_ext)
                #
                if not verbose: continue
                domain_residue = [(1+np.where(self.segment_index == segment_i[l])[0]).tolist() for l in sorted(xi)]
                domain = []
                for rr in domain_residue:
                    if len(domain) == 0:
                        domain.append(rr)
                    else:
                        if domain[-1][-1]+1 == rr[0]:
                            domain[-1].extend(rr)
                        else:
                            domain.append(rr)
                print (":"+','.join(['%d-%d'%(d[0],d[-1]) for d in domain]))
            return domain_s
        #
        domain_out = []
        #
        sse_map = copy.deepcopy(self.map0)
        mask = (self.segment_index >= 0)
        prob_sse = prob[np.ix_(mask,mask)]
        #
        seg = self.segment_index[mask]
        np.add.at(sse_map, np.ix_(seg,seg), prob_sse)
        #
        has_contact = (sse_map > PARAM_PROB_CUTOFF)
        n_graphs, labels = csgraph.connected_components(has_contact)
        #
        segment_no = np.arange(self.n_segment)
        for i in range(n_graphs):
            sub = (labels==i)
            if np.where(sub)[0].shape[0] < PARAM_DOMAIN_MIN_SEG:
                continue
            #
            segment_i = segment_no[sub]
            graph = nx.Graph(has_contact[np.ix_(sub, sub)])
            community = nx.algorithms.community.girvan_newman(graph)
            #
            Xs = [np.arange(len(segment_i)).astype(int)]
            domain0 = get_domain(self, Xs, segment_i, verbose=verbose)
            domain_size = [np.where(d==True)[0].shape[0] for d in domain0]
            if min(domain_size) < PARAM_DOMAIN_MIN_RES:
                continue
            #
            for X in community:
                nX = [len(xi) for xi in X]
                if min(nX) < PARAM_DOMAIN_MIN_SEG:
                    break
                #
                Xi = []
                for xi in X:
                    Xi.append(segment_i[np.array(list(xi), dtype=int)])
                nXi = len(Xi)
                n_contact_between_sub = np.eye(nXi, dtype=int)
                for a in range(nXi-1):
                    for b in range(a+1, nXi):
                        map_ab = has_contact[np.ix_(Xi[a], Xi[b])]
                        n_contact_between_sub[a,b] = map_ab.astype(int).sum()
                n_contact_between_sub += n_contact_between_sub.T
                if n_contact_between_sub.max() > PARAM_DOMAIN_CONTACT_MIN:
                    break

                domain = get_domain(self, X, segment_i, verbose=verbose)
                domain_size = [min((np.where(d==True)[0].shape[0]), np.where(d!=True)[0].shape[0]) for d in domain]
                if min(domain_size) < PARAM_DOMAIN_MIN_RES:
                    break
                else:
                    domain0 = domain
                break

            domain_out.extend(domain0)
        return domain_out

def run(ss, contact, verbose=False):
    segment = SegmentMap(ss)
    domain_new = segment.apply_map(contact, verbose=verbose)
    return domain_new

def read_psipred(psipred_fn):
    ss = []
    with open(psipred_fn) as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            line = line.strip()
            if len(line) == 0:
                continue
            ss.append(line.split()[2])
    ss = ''.join(ss)
    return ss

def get_soft_domain_boundary(mask0, soft=PARAM_DOMAIN_BOUNDARY*2):
    mask = np.ones(mask0.shape[0], dtype=float)
    #
    weights = np.arange(1, 1+soft, dtype=float) / (1.+soft)
    #
    unmask = np.array(np.where(~mask0)[0], dtype=int)
    mask[~mask0] = 0.0
    for i,x in enumerate(mask0):
        if not x: continue
        d = np.min(np.abs(unmask-i))
        if d > soft: continue
        mask[i] = weights[d-1]
    return mask

def get_contact(npz_fn):
    npz = np.load(npz_fn.short())
    prob = npz['dist'][:,:,1:13].sum(axis=-1)
    prob[prob < 0.15] = 0.0
    return prob

def main():
    ss = read_psipred(path.Path(sys.argv[1]).short())
    npz = get_contact(path.Path(sys.argv[2]))
    run(ss, npz, verbose=True)

if __name__ == '__main__':
    main()
