#!/bin/bash

RWHOME=/apps/rwplus
EXEC=./calRWplus

n=$#

if [ $n -eq 0 ]
then
    echo "USAGE: calc_rwplus.sh [PDBs]"
    exit -1
fi


pwd=$(pwd)
for pdb_fn in $@;
do
    abs_fn=$(readlink -f $pdb_fn)
    cd $RWHOME
    score=$($EXEC $abs_fn | cut -d= -f2 | awk '{print $1}')
    printf "%12.4f  %s\n" $score $pdb_fn
    cd $pwd
done
