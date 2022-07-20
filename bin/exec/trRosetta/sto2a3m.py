#!/usr/bin/env python

import sys
import path


def sto2a3m(sto):
    name_s = []
    seq_s = []
    with sto.open() as fp:
        n_chunk = -1
        for line in fp:
            if line.startswith("#"):
                continue
            if line.startswith("//"):
                break
            if line.strip() == "":
                i = -1
                continue
            i += 1
            if i == 0:
                n_chunk += 1
            #
            x = line.strip().split()
            name = x[0]
            seq = x[1]
            if n_chunk == 0:
                name_s.append(name)
                seq_s.append("")
            seq_s[i] += seq
    #
    for i in range(len(seq_s)):
        seq_s[i] = seq_s[i].replace(".", "")
    #
    a3m = []
    seqq = seq_s[0]
    for seqt in seq_s:
        a3m_seq = []
        for q, t in zip(seqq, seqt):
            if q != "-":
                a3m_seq.append(t)
            elif t != "-":
                a3m_seq.append(t.lower())
        a3m.append("".join(a3m_seq))

    wrt = []
    for i, name in enumerate(name_s):
        wrt.append(">%s\n" % name)
        wrt.append(a3m[i] + "\n")

    return wrt


def main():
    sto_fn = path.Path(sys.argv[1])
    a3m = sto2a3m(sto_fn)
    a3m_fn = path.Path("%s.a3m" % (sto_fn.prefix()))
    with a3m_fn.open("wt") as fout:
        fout.writelines(a3m)


if __name__ == "__main__":
    main()
