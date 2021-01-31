#!/bin/bash
FILE=`readlink -e $0`
rootdir=`dirname $FILE`
cd $rootdir/prediction
for species in `cat $rootdir/input/list`;do
    echo "predict GO for species $species"
    $rootdir/bin/predict_nw.py $rootdir/input/target.$species.fasta _${species}_go.txt
done
