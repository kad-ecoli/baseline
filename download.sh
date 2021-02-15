#!/bin/bash
FILE=`readlink -e $0`
rootdir=`dirname $FILE`

cd $rootdir/download

echo download CAFA3 data
wget https://www.biofunctionprediction.org/cafa-targets/CAFA3_targets.tgz -O CAFA3_targets.tgz
wget https://www.biofunctionprediction.org/cafa-targets/CAFA3_training_data.tgz -O CAFA3_training_data.tgz
wget https://ndownloader.figshare.com/files/17519846 -O supplementary_data.tar.gz
wget https://www.biofunctionprediction.org/annotations/gene_ontology_edit.obo.2016-06-01.gz -O gene_ontology_edit.obo.2016-06-01.gz
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+-x64-linux.tar.gz -O ncbi-blast-2.4.0+-x64-linux.tar.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/goa_uniprot_all.gaf.157.gz -O goa_uniprot_all.gaf.157.gz

echo decompress CAFA3 data
tar -xvf CAFA3_targets.tgz
tar -xvf CAFA3_training_data.tgz
tar -xvf supplementary_data.tar.gz
tar -xvf supplementary_data/cafa3/benchmark20171115.tar
tar -xvf ncbi-blast-2.4.0+-x64-linux.tar.gz

echo process CAFA3 data
cp $rootdir/download/ncbi-blast-2.4.0+/bin/makeblastdb $rootdir/bin
cp $rootdir/download/ncbi-blast-2.4.0+/bin/blastp      $rootdir/bin
cp $rootdir/download/ncbi-blast-2.4.0+/bin/blastdbcmd  $rootdir/bin
cat $rootdir/download/benchmark20171115/lists/*o_all_type1.txt |sort |uniq > $rootdir/download/benchmark20171115/lists/all_type1.txt
cat $rootdir/download/benchmark20171115/lists/*o_all_type2.txt |sort |uniq > $rootdir/download/benchmark20171115/lists/all_type2.txt
cat $rootdir/download/benchmark20171115/lists/all_type1.txt $rootdir/download/benchmark20171115/lists/all_type2.txt > $rootdir/download/benchmark20171115/lists/all_type.txt 

cd "$rootdir/download/Target files"
ls target*fasta > list
for prefix in `cat list|sed 's/.fasta//g'`;do
    $rootdir/bin/makeblastdb -in $prefix.fasta -dbtype prot -parse_seqids
    $rootdir/bin/blastdbcmd  -db $prefix.fasta -dbtype prot -entry_batch $rootdir/download/benchmark20171115/lists/all_type.txt -out $rootdir/input/$prefix.fasta
    if [ -f "$rootdir/input/$prefix.fasta" ] && [ ! -s "$rootdir/input/$prefix.fasta" ];then
	rm "$rootdir/input/$prefix.fasta"
    fi
    rm $prefix.fasta.p*
done
ls $rootdir/input/target.*.fasta|cut -f2 -d. > $rootdir/input/list

cd $rootdir/data
zcat $rootdir/download/gene_ontology_edit.obo.2016-06-01.gz > $rootdir/data/go-basic.obo
cp $rootdir/download/CAFA3_training_data/uniprot_sprot_exp.fasta $rootdir/data/uniprot_sprot_exp.fasta
$rootdir/bin/makeblastdb -in uniprot_sprot_exp.fasta -dbtype prot -parse_seqids
$rootdir/bin/propagate_training_terms.py go-basic.obo $rootdir/download/CAFA3_training_data/uniprot_sprot_exp.txt
zcat $rootdir/download/goa_uniprot_all.gaf.157.gz|grep -P `cat $rootdir/input/list |sed 's/^/taxon:/g'|paste -sd'|'|sed 's/^/(/g'|sed 's/$/\b)/g'`|grep -P "^UniProtKB" |grep -vP "\tND\t"|grep -vP "\tNOT\b"  > $rootdir/data/goa_uniprot_all.gaf
$rootdir/bin/cull_IEA.py go-basic.obo $rootdir/input/target.map goa_uniprot_all.gaf goa_uniprot_all.clean.gaf goa_uniprot_all.is_a
mv goa_uniprot_all.clean.gaf goa_uniprot_all.gaf 

cd $rootdir/groundtruth
$rootdir/bin/propagate_groundtruth_terms.py $rootdir/data/go-basic.obo F $rootdir/download/benchmark20171115/groundtruth/leafonly_MFO.txt leafonly_MFO.is_a 
$rootdir/bin/propagate_groundtruth_terms.py $rootdir/data/go-basic.obo P $rootdir/download/benchmark20171115/groundtruth/leafonly_BPO.txt leafonly_BPO.is_a 
$rootdir/bin/propagate_groundtruth_terms.py $rootdir/data/go-basic.obo C $rootdir/download/benchmark20171115/groundtruth/leafonly_CCO.txt leafonly_CCO.is_a 
cp $rootdir/download/benchmark20171115/lists/all_type1.txt .
cp $rootdir/download/benchmark20171115/lists/all_type2.txt .
