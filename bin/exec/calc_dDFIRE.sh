#!/bin/bash

export DATADIR=$HOME/apps/dfire/
EXEC=$HOME/apps/dfire/dDFIRE

n=$#

if [ $n -eq 0 ]
then
    echo "USAGE: calc_dDFIRE.sh [PDBs]"
    exit -1
fi


pwd=$(pwd)
for pdb_fn in $@;
do
    abs_fn=$(readlink -f $pdb_fn)
    score=$($EXEC $abs_fn | cut -d: -f2 | awk '{print $1, $2}')
    printf "%8.3f  %8.3f  %s\n" $score $pdb_fn
done
