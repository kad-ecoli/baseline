#!/usr/bin/env python
docstring='''predict_blastbitscore.py target.9606.fasta blastbitscore_9606_go.txt
    use blast bitscore to predict GO terms for input file target.9606.fasta,
    output result to blastbitscore_9606_go.txt
'''
import sys
import os
from os.path import dirname, basename, abspath
from subprocess import Popen,PIPE
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

def run_self_blast(infile,outfile):
    cmd="cp %s %s.fasta; %s/makeblastdb -in %s.fasta -dbtype prot -parse_seqids"%(infile,outfile,bindir,outfile)
    os.system(cmd)
    cmd="%s/blastp -db %s.fasta -outfmt '6 qacc sacc bitscore' -query %s"%(bindir,outfile,infile)
    p=Popen(cmd,shell=True,stdout=PIPE)
    stdout,stderr=p.communicate()
    selfscore_dict=dict()
    for line in stdout.splitlines():
        qacc,sacc,bitscore=line.split('\t')
        if qacc!=sacc:
            continue
        selfscore_dict[qacc]=float(bitscore)
    cmd="rm %s.fasta*"%(outfile)
    os.system(cmd)
    return selfscore_dict

def run_blast(infile,selfscore_dict):
    cmd="%s/blastp -db %s/uniprot_sprot_exp.fasta -outfmt '6 qacc sacc bitscore' -query %s"%(bindir,datdir,infile)
    p=Popen(cmd,shell=True,stdout=PIPE)
    stdout,stderr=p.communicate()
    blast_dict=dict()
    for line in stdout.splitlines():
        qacc,sacc,bitscore=line.split('\t')
        if not qacc in selfscore_dict:
            print("ERROR! missing self score for "+qacc)
            continue
        bitscore=float(bitscore)
        if not qacc in blast_dict:
            blast_dict[qacc]=[]
        blast_dict[qacc].append([sacc,min(1,bitscore/selfscore_dict[qacc])])
    return blast_dict

def write_output(blast_dict,annotation_dict,outfile):
    print("writing "+outfile)
    txt ="AUTHOR blastbitscore\n"
    txt+="MODEL 1\n"
    txt+="KEYWORDS sequence alignment.\n"
    for target in sorted(blast_dict.keys()):
        for Aspect in "FPC":
            predict_dict=dict()
            denominator=0
            for sacc,score in blast_dict[target]:
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
    outfile=sys.argv[2]
    selfscore_dict =run_self_blast(infile,outfile)
    annotation_dict=read_annotation()
    blast_dict     =run_blast(infile,selfscore_dict)
    write_output(blast_dict,annotation_dict,outfile)
