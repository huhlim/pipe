#!/usr/bin/env python

import os
import sys
import glob
import sqlite3
try:
    from rich.progress import track
except:
    def track(X, description='Running...'):
        for x in X:
            yield x

PROTEOME_DIR='/mnt/ramdisk/proteome.dir'
#PROTEOME_DIR='/green/s2/huhlim/db/uniref/proteome.dir'

def read_a3m(a3m_fn):
    a3m = [[], []]
    with open(a3m_fn) as fp:
        for line in fp:
            if line.startswith(">"):
                name = line.strip().split()[0][1:]
                if len(a3m[0]) == 0:
                    uniref_id = None
                else:
                    uniref_id = name.split("_")[1].split("/")[0]
                a3m[0].append([name, uniref_id, None])
            else:
                a3m[1].append(line)
    return a3m

def write_a3m(a3m):
    wrt = []
    for i in range(len(a3m[0])):
        if a3m[0][i][1] is None:
            wrt.append(">%s\n"%a3m[0][i][0])
            wrt.append("%s"%a3m[1][i])
        if a3m[0][i][2] is None:
            continue
        wrt.append(">%s %s\n"%(a3m[0][i][0], a3m[0][i][2]))
        wrt.append("%s"%a3m[1][i])
    return wrt

def get_proteome(a3m):
    db_s = [x.split("/")[-1][:-3] for x in glob.glob("%s/*.db"%(PROTEOME_DIR))]
    #
    key_s = {}
    for i,entry in enumerate(a3m[0]):
        uniref_id = entry[1]
        if uniref_id is None:
            continue
        for key_len in [5,4,3]:
            key = uniref_id[:key_len]
            if key in db_s:
                break
        if key not in key_s:
            key_s[key] = []
        key_s[key].append((i,uniref_id))
    #
    for key in track(key_s):
        db_fn = '%s/%s.db'%(PROTEOME_DIR, key)
        if not os.path.exists(db_fn):
            continue
        #
        conn = sqlite3.connect(db_fn)
        cursor = conn.cursor()
        for index, uniref_id in key_s[key]:
            cursor.execute("SELECT (proteome) FROM mapping where id='%s' LIMIT 1;"%(uniref_id))
            proteome = cursor.fetchone()
            if proteome is not None:
                a3m[0][index][2] = proteome[0].replace(" ","")
        conn.close()
    #
    return a3m

def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: %s [in_A3M] [out_A3M]\n"%__file__)
    #
    in_a3m = sys.argv[1]
    out_a3m = sys.argv[2]
    a3m = read_a3m(in_a3m)
    #
    get_proteome(a3m)
    #
    with open(out_a3m, 'wt') as fout:
        fout.writelines(write_a3m(a3m))

if __name__ == '__main__':
    main()
