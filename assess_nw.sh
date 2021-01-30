#!/bin/bash
FILE=`readlink -e $0`
rootdir=`dirname $FILE`
cd $rootdir/CAFA_assessment_tool

method_list="nwlocalID_1
nwalnscore_1
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
    outfile="$rootdir/CAFA_assessment_tool/${method}_all.txt"
    head -3 $infile > $outfile
    for species in `cat $rootdir/input/list`;do
        infile="$rootdir/prediction/${method}_${species}_go.txt"
    	head -n-1 $infile |tail -n+4 >> $outfile
    done
    tail -1 $infile >> $outfile
    if [ ! -s "$outfile" ];then
        continue
    fi
    outfile=`basename $outfile`
    configfile="$rootdir/CAFA_assessment_tool/${method}_all.yaml"
    echo "---
assess:
    file: $outfile
    obo: ./precrec/go_cafa3.obo
    benchmark: ./precrec/benchmark/CAFA3_benchmarks/
    results: ./results/"   > $configfile
    python3 $rootdir/CAFA_assessment_tool/assess_main.py $configfile
done
