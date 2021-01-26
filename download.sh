#!/bin/bash
FILE=`readlink -e $0`
rootdir=`dirname $FILE`
cd $rootdir

echo download CAFA3 data

wget https://www.biofunctionprediction.org/cafa-targets/CAFA3_targets.tgz -O $rootdir/download/CAFA3_targets.tgz
wget https://www.biofunctionprediction.org/cafa-targets/CAFA3_training_data.tgz -O $rootdir/download/CAFA3_training_data.tgz
wget https://ndownloader.figshare.com/files/17519846 -O $rootdir/download/supplementary_data.tar.gz
wget https://www.biofunctionprediction.org/annotations/gene_ontology_edit.obo.2016-06-01.gz -O $rootdir/download/gene_ontology_edit.obo.2016-06-01.gz

echo decompress CAFA3 data

cd $rootdir/download
tar -xvf CAFA3_targets.tgz
tar -xvf CAFA3_training_data.tgz
tar -xvf supplementary_data.tar.gz
tar -xvf supplementary_data/cafa3/benchmark20171115.tar
