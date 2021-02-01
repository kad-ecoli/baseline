#!/bin/bash
FILE=`readlink -e $0`
rootdir=`dirname $FILE`
cd $rootdir/prediction
for species in `cat $rootdir/input/list`;do
    echo "predict GO for species $species"
    $rootdir/bin/predict_iea.py $rootdir/input/target.$species.fasta iea_1_${species}_go.txt
    $rootdir/bin/predict_naive.py $rootdir/input/target.$species.fasta naive_1_${species}_go.txt
    $rootdir/bin/predict_blast.py $rootdir/input/target.$species.fasta _${species}_go.txt
    $rootdir/bin/predict_blastbitscore.py $rootdir/input/target.$species.fasta blastbitscore_1_${species}_go.txt
done
