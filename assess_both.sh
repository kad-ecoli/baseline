#!/bin/bash
FILE=`readlink -e $0`
rootdir=`dirname $FILE`
cd $rootdir/results

method_list="naive_1
iea_1
blastlocalID_1
blastglobalID_1
blastglobalID_2
blastglobalID_3
blastevalue_1
blastrank_1
blastfreq_1
blastmetago_1
blastnetgo_1
blastbitscore_1
blastbitscore_2
blastbitscore_3

nwlocalID_1
nwalnscore_1
nwalnscore_2
nwalnscore_3
nwglobalID_1
nwglobalID_2
nwglobalID_3
nwrank_1
nwfreq_1
nwmetago_1
nwnetgo_1"

for method in $method_list;do
    echo assessing $method
    species=`head -1 $rootdir/input/list`
    infile="$rootdir/prediction/${method}_${species}_go.txt"
    outfile="$rootdir/results/${method}_all.txt"
    resultfile="$rootdir/results/${method}_all_results.txt"
    head -3 $infile > $outfile
    for species in `cat $rootdir/input/list`;do
        infile="$rootdir/prediction/${method}_${species}_go.txt"
    	head -n-1 $infile |tail -n+4 >> $outfile
    done
    tail -1 $infile >> $outfile
    if [ ! -s "$outfile" ];then
        continue
    fi
    $rootdir/bin/assess_result.py $outfile $resultfile
done
./plot.py
./plot_nw.py
