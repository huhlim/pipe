#!/usr/bin/env python

import os
import sys


def read_a3m(fn):
    name = []
    sequence = []
    with open(fn) as fp:
        for line in fp:
            if line.startswith(">"):
                name.append(line.strip())
                seq = []
                sequence.append(seq)
            else:
                seq.append(line.strip())
    for i in range(len(sequence)):
        sequence[i] = "".join(sequence[i])
    return name, sequence


def filter_a3m(name_s, sequence_s):
    out_s = [[], []]
    for name, sequence in zip(name_s, sequence_s):
        status = True
        for i, out in enumerate(out_s[0]):
            if out == name and out_s[1][i] == sequence:
                status = False
                break
        if status:
            out_s[0].append(name)
            out_s[1].append(sequence)
    return out_s


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: %s [in A3M] [out A3M]" % (__file__))
    #
    in_A3M = sys.argv[1]
    out_A3M = sys.argv[2]
    #
    name_s, sequence_s = read_a3m(in_A3M)
    name_s, sequence_s = filter_a3m(name_s, sequence_s)
    #
    with open(out_A3M, "wt") as fout:
        for name, seq in zip(name_s, sequence_s):
            fout.write("%s\n" % name)
            fout.write("%s\n" % seq)


if __name__ == "__main__":
    main()
