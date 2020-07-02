#!/usr/bin/env python

import sys

def convert_to_res(Xs):
    res_s = []
    res = [] ; res_s.append(res)
    for i,x in enumerate(Xs):
        if x == 'X':
            if len(res) > 1 and i-1 not in res:
                res = [] ; res_s.append(res)
            res.append(i)

    out = []
    for res in res_s:
        out.append("%d-%d"%(res[0]+1, res[-1]+1))
    return ','.join(out)

def main():
    with open(sys.argv[1]) as fp:
        for line in fp:
            domain_name, Xs = line.strip().split()
            sys.stdout.write("%-14s  %s\n"%(domain_name, convert_to_res(Xs)))

if __name__ == '__main__':
    main()
