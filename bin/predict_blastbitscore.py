#!/usr/bin/env python
docstring='''predict_blastbitscore.py target.9606.fasta _9606_go.txt
    use blast bitscore to predict GO terms for input file target.9606.fasta,
output:
    blastbitscore_1_9606_go.txt - bitscore / selfscore(query)
    blastbitscore_2_9606_go.txt - bitscore / selfscore(subject)
    blastbitscore_3_9606_go.txt - bitscore / max{selfscore(query),selfscore(subject)}
'''
import sys
import os
from os.path import dirname, basename, abspath
from subprocess import Popen,PIPE
from blastSelfScore import run_self_blast_all
bindir=dirname(abspath(__file__))
rootdir=dirname(bindir)
datdir=rootdir+"/data"

def read_annotation():
    annotation_dict=dict()
    for Aspect in "FPC":
        annotation_dict[Aspect]=dict()
        filename=datdir+"/uniprot_sprot_exp."+Aspect
        fp=open(filename,'r')
        for line in fp.read().splitlines():
            target,GOterm_list=line.split('\t')
            annotation_dict[Aspect][target]=GOterm_list.split(',')
        fp.close()
    return annotation_dict

def read_db_selfscore():
    dbscore_dict=dict()
    filename=datdir+"/uniprot_sprot_exp.fasta.SelfScore_blast"
    fp=open(filename,'r')
    for line in fp.read().splitlines():
        target,bitscore=line.split('\t')
        dbscore_dict[target]=float(bitscore)
    fp.close()
    return dbscore_dict

def run_blast(infile,selfscore_dict,dbscore_dict):
    cmd="%s/blastp -db %s/uniprot_sprot_exp.fasta -outfmt '6 qacc sacc bitscore' -query %s"%(bindir,datdir,infile)
    p=Popen(cmd,shell=True,stdout=PIPE)
    stdout,stderr=p.communicate()
    blast_dict=dict()
    for line in stdout.splitlines():
        qacc,sacc,bitscore=line.split('\t')
        bitscore=float(bitscore)
        if not qacc in blast_dict:
            blast_dict[qacc]=[]
        blast_dict[qacc].append([sacc,
            min(1,bitscore/selfscore_dict[qacc]),
            min(1,bitscore/dbscore_dict[sacc]),
            min(1,bitscore/max(selfscore_dict[qacc],dbscore_dict[sacc])),
            ])
    return blast_dict

def write_output(blast_dict,annotation_dict,suffix):
    for m in [1,2,3]:
        outfile="blastbitscore_%d%s"%(m,suffix)
        print("writing "+outfile)
        txt ="AUTHOR blastbitscore\n"
        txt+="MODEL %d\n"%m
        txt+="KEYWORDS sequence alignment.\n"
        for target in sorted(blast_dict.keys()):
            for Aspect in "FPC":
                predict_dict=dict()
                denominator=0
                for line in blast_dict[target]:
                    sacc=line[0]
                    score=line[m]
                    if not sacc in annotation_dict[Aspect]:
                        continue
                    GOterm_list=annotation_dict[Aspect][sacc]
                    for GOterm in GOterm_list:
                        if not GOterm in predict_dict or \
                            score>predict_dict[GOterm]:
                            predict_dict[GOterm]=score
                for cscore,GOterm in sorted([(predict_dict[GOterm],
                    GOterm) for GOterm in predict_dict],reverse=True):
                    cscore="%.2f"%cscore
                    if cscore=="0.00":
                        break
                    txt+="%s\t%s\t%s\n"%(target,GOterm,cscore)
        txt+="END\n"
        fp=open(outfile,'w')
        fp.write(txt)
        fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=3:
        sys.stderr.write(docstring)
        exit()

    infile=sys.argv[1]
    suffix=sys.argv[2]
    target_list,selfscore_dict=run_self_blast_all(infile,"blastbitscore"+suffix)
    annotation_dict=read_annotation()
    dbscore_dict   =read_db_selfscore()
    blast_dict     =run_blast(infile,selfscore_dict,dbscore_dict)
    write_output(blast_dict,annotation_dict,suffix)
