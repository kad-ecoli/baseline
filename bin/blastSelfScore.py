#!/usr/bin/env python
import os
docstring='''blastSelfScore.py target.9606.fasta target.9606.fasta.SelfScore_blast
    compute the blastp bit-score for self-alignment
output:
    target.9606.fasta.SelfScore_blast
'''
import sys
from os.path import dirname, basename, abspath
from subprocess import Popen,PIPE
from math import exp
bindir=dirname(abspath(__file__))
rootdir=dirname(bindir)
batchsize=10

def read_sequence_as_dict(infile):
    target_list=[]
    sequence_dict=dict()
    fp=open(infile,'r')
    for block in ('\n'+fp.read()).split('\n>'):
        lines=block.strip().splitlines()
        if len(lines)<2:
            continue
        target=lines[0].split()[0]
        sequence=''.join(lines[1:])
        target_list.append(target)
        sequence_dict[target]=sequence
    fp.close()
    return target_list,sequence_dict

def run_blast_cmd(blastp,outfile,selfscore_dict):
    cmd="%s/makeblastdb -in %s.fasta -dbtype prot -parse_seqids"%(bindir,outfile)
    os.system(cmd)
    cmd="%s -db %s.fasta -outfmt '6 qacc sacc bitscore' -query %s.fasta"%(blastp,outfile,outfile)
    p=Popen(cmd,shell=True,stdout=PIPE)
    stdout,stderr=p.communicate()
    for line in stdout.splitlines():
        qacc,sacc,bitscore=line.split('\t')
        if qacc!=sacc:
            continue
        selfscore_dict[qacc]=float(bitscore)
    cmd="rm %s.fasta*"%(outfile)
    os.system(cmd)
    return selfscore_dict
    
def make_unscored_fasta(target_list,sequence_dict,selfscore_dict,outfile):
    txt=''
    for target in target_list:
        if not target in selfscore_dict:
            sys.stderr.write("scoring %s\n"%target)
            txt+=">%s\n%s\n"%(target,sequence_dict[target])
    if len(txt)==0:
        return False

    fp=open(outfile+".fasta",'w')
    fp.write(txt)
    fp.close()
    return True

def run_self_blast(infile,outfile,target_list,sequence_dict):
    selfscore_dict=dict()
    make_unscored_fasta(target_list,sequence_dict,selfscore_dict,outfile)

    blastp="%s/blastp -max_target_seqs 1"%bindir
    selfscore_dict=run_blast_cmd(blastp,outfile,selfscore_dict)
    if not make_unscored_fasta(target_list,sequence_dict,selfscore_dict,outfile):
        return selfscore_dict
    
    blastp="%s/blastp"%bindir
    selfscore_dict=run_blast_cmd(blastp,outfile,selfscore_dict)
    if not make_unscored_fasta(target_list,sequence_dict,selfscore_dict,outfile):
        return selfscore_dict
    
    blastp="%s/blastp -task blastp-short"%bindir
    selfscore_dict=run_blast_cmd(blastp,outfile,selfscore_dict)
    return selfscore_dict

def fallback_SelfScore(infile):
    selfscore_all_dict=dict()
    cmd="%s/SelfScore %s"%(bindir,infile)
    p=Popen(cmd,shell=True,stdout=PIPE)
    stdout,stderr=p.communicate()
    for line in stdout.splitlines():
        if line.startswith('#'):
            continue
        qacc,length,score,bitscore=line.split('\t')
        selfscore_all_dict[qacc]=float(bitscore)
    return selfscore_all_dict

def write_output(target_list,selfscore_dict,outfile):
    txt=''
    for target in target_list:
        if target in selfscore_dict:
            txt+="%s\t%s\n"%(target,selfscore_dict[target])
    fp=open(outfile,'w')
    fp.write(txt)
    fp.close()
    return

def run_self_blast_all(infile,outfile):
    target_all_list,sequence_dict=read_sequence_as_dict(infile)
    selfscore_all_dict=fallback_SelfScore(infile)
    for i in range(0,len(target_all_list),batchsize):
        target_list=target_all_list[i:i+batchsize]
        selfscore_dict=run_self_blast(infile,outfile,target_list,sequence_dict)
        for target in selfscore_dict:
            selfscore_all_dict[target]=selfscore_dict[target]
    return target_all_list,selfscore_all_dict

if __name__=="__main__":
    if len(sys.argv)!=3:
        sys.stderr.write(docstring)
        exit()

    infile=sys.argv[1]
    outfile=sys.argv[2]
    target_all_list,selfscore_all_dict=run_self_blast_all(infile,outfile)
    write_output(target_all_list,selfscore_all_dict,outfile)
