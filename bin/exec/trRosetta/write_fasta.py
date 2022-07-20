#!/usr/bin/env python

import sys
import path


def main():
    if len(sys.argv) == 1:
        sys.exit("usage: %s [FA]" % __file__)
    #
    fa_fn = path.Path(sys.argv[1])
    title = fa_fn.name()
    seq = []
    with fa_fn.open() as fp:
        for line in fp:
            if not line.startswith(">"):
                seq.append(line.strip())
    seq = "".join(seq)
    #
    domain_s = {}
    with open("domain_s") as fp:
        for line in fp:
            name, mask_str = line.strip().split()
            mask = [X == "X" for X in mask_str]
            domain_s[name] = mask
    #
    for domain_name in domain_s:
        fn = path.Path("%s/%s.fa" % (domain_name, title))
        if fn.status():
            continue
        domain_seq = "".join([s for s, m in zip(seq, domain_s[domain_name]) if m])
        #
        with fn.open("wt") as fout:
            fout.write(">%s\n" % title)
            fout.write(domain_seq)


if __name__ == "__main__":
    main()
