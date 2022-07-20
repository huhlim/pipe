#!/usr/bin/env python

import os
import sys

try:
    from rich.progress import track
except:

    def track(X, description="Running..."):
        for x in X:
            yield x


SPACERq = "G" * 20
SPACERt = "-" * (len(SPACERq))


class Sequence(object):
    def __init__(self, header, is_query=False):
        self.has_pair = False
        self.is_query = is_query
        self.header = header.split(maxsplit=1)
        if self.is_query:
            self.prot_s = set([])
        else:
            self.prot_s = set(self.header[1].replace(" ", "").split(","))
        self._seq = []

    def __len__(self):
        return len(self.seq)

    def __eq__(self, othr):
        return len(self.prot_s.intersection(othr.prot_s)) > 0

    def append(self, line):
        self._seq.append(line)

    @property
    def seq(self):
        return "".join(self._seq)

    def __repr__(self):
        return self.seq

    @classmethod
    def create_empty_sequence(cls, seq):
        s = cls("empty", is_query=True)
        s._seq = "-" * (len(seq))
        return s


class PairedSequence(object):
    def __init__(self, i, seq_1, seq_2):
        self.header = ">pair_%d" % (i)
        if not (seq_1.is_query and seq_2.is_query):
            prot = list(seq_1.prot_s.intersection(seq_2.prot_s))
            prot = ",".join(prot)
            self.header = "%s %s %s %s" % (
                self.header,
                prot,
                seq_1.header[0].replace(">", ""),
                seq_2.header[0].replace(">", ""),
            )
            self.seq = "%s%s%s" % (seq_1.seq, SPACERt, seq_2.seq)
        else:
            self.seq = "%s%s%s" % (seq_1.seq, SPACERq, seq_2.seq)
        seq_1.has_pair = True
        seq_2.has_pair = True

    def write(self):
        return self.__repr__()

    def __repr__(self):
        return "%s\n%s\n" % (self.header, self.seq)


def read_a3m(fn):
    a3m = []
    with open(fn) as fp:
        for line in fp:
            if line.startswith(">"):
                seq = Sequence(line.strip(), is_query=(len(a3m) == 0))
                a3m.append(seq)
            else:
                seq.append(line.strip())
    return a3m


def pair_a3m(a3m_1, a3m_2, paired_only):
    pair_s = []
    pair_s.append(PairedSequence(0, a3m_1[0], a3m_2[0]))
    index = 0
    for seq_1 in track(a3m_1[1:]):
        for seq_2 in a3m_2[1:]:
            if seq_1 == seq_2:
                index += 1
                pair_s.append(PairedSequence(index, seq_1, seq_2))
    #
    if paired_only:
        return pair_s
    empty_1 = Sequence.create_empty_sequence(a3m_1[0])
    empty_2 = Sequence.create_empty_sequence(a3m_2[0])
    #
    for seq_1 in a3m_1[1:]:
        if seq_1.has_pair:
            continue
        index += 1
        pair_s.append(PairedSequence(index, seq_1, empty_2))
    for seq_2 in a3m_2[1:]:
        if seq_2.has_pair:
            continue
        index += 1
        pair_s.append(PairedSequence(index, empty_1, seq_2))

    return pair_s


def main():
    if len(sys.argv) < 4:
        sys.exit("USAGE: %s [A3M_1] [A3M_2] [out_A3M]\n" % __file__)
    #
    a3m_1_fn = sys.argv[1]
    a3m_2_fn = sys.argv[2]
    out_a3m = sys.argv[3]
    if len(sys.argv) > 4 and sys.argv[4] == "--pair":
        paired_only = True
    else:
        paired_only = False
    #
    a3m_1 = read_a3m(a3m_1_fn)
    a3m_2 = read_a3m(a3m_2_fn)
    #
    paired = pair_a3m(a3m_1, a3m_2, paired_only)
    with open(out_a3m, "wt") as fout:
        for pair in paired:
            fout.write(pair.write())


if __name__ == "__main__":
    main()
