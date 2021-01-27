#!/bin/bash
FILE=`readlink -e $0`
rootdir=`dirname $FILE`

for species in `cat $rootdir/input/list`;do
    echo "predict GO for species $species"
    $rootdir/bin/predict_naive.py $rootdir/input/target.$species.fasta $rootdir/prediction/naive_1_${species}_go.txt
done
