#!/usr/bin/env python
docstring='''predict_blastbest.py target.9606.fasta blastbest_1_9606_go.txt
    combine blast and naive run to predict GO terms for input file
    target.9606.fasta, output prediction to blastbest_1_9606_go.txt
'''
import sys
from os.path import dirname, basename, abspath
from subprocess import Popen,PIPE
from math import exp
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

def read_naive_prob():
    naive_GOterm_dict=dict()
    for Aspect in "FPC":
        naive_GOterm_dict[Aspect]=[]
        filename=datdir+"/naive."+Aspect
        fp=open(filename,'r')
        for line in fp.read().splitlines():
            GOterm,cscore=line.split('\t')[:2]
            cscore=float(cscore)
            if "%.2f"%cscore=="0.00":
                break
            naive_GOterm_dict[Aspect].append((GOterm,cscore))
        fp.close()
    return naive_GOterm_dict

def read_fasta_as_list(infile):
    target_list=[]
    fp=open(infile,'r')
    for line in fp.read().splitlines():
        if line.startswith('>'):
            target_list.append(line[1:].split()[0])
    fp.close()
    return target_list

def run_blast(infile):
    cmd="%s/blastp -db %s/uniprot_sprot_exp.fasta -outfmt '6 qacc qlen sacc slen bitscore nident' -query %s"%(bindir,datdir,infile)
    p=Popen(cmd,shell=True,stdout=PIPE)
    stdout,stderr=p.communicate()
    blast_dict=dict()
    for line in stdout.splitlines():
        qacc,qlen,sacc,slen,bitscore,nident=line.split('\t')
        qlen=float(qlen)
        slen=float(slen)
        bitscore=float(bitscore)
        nident=float(nident)
        if not qacc in blast_dict:
            blast_dict[qacc]=[]
        #blast_dict[qacc].append([sacc, 2*nident/(qlen+slen), bitscore])
        blast_dict[qacc].append([sacc, nident*(1/qlen + 1/slen)/2, bitscore])
    return blast_dict

def write_output(target_list,blast_dict,annotation_dict,naive_GOterm_dict,outfile):
    txt ="AUTHOR blastbest\n"
    txt+="MODEL 1\n"
    txt+="KEYWORDS sequence alignment.\n"
    for target in target_list:
        if not target in blast_dict:
            for Aspect in "FPC":
                for GOterm,cscore in naive_GOterm_dict[Aspect]:
                    txt+="%s\t%s\t%.2f\n"%(target,GOterm,cscore)
            continue
        for Aspect in "FPC":
            predict_dict=dict()
            denominator=0.
            total_score=1.
            #template_num=0
            for template,globalID,score in blast_dict[target]:
                if not template in annotation_dict[Aspect]:
                    continue
                #template_num+=1
                total_score*=(1-globalID)
                denominator+=score
                GOterm_list=annotation_dict[Aspect][template]
                for GOterm in GOterm_list:
                    if not GOterm in predict_dict:
                        predict_dict[GOterm]=0.
                    predict_dict[GOterm]+=score
            total_score=1-total_score
            #print(target,Aspect,total_score,template_num)
            for GOterm in predict_dict:
                #predict_dict[GOterm]*=total_score/denominator
                predict_dict[GOterm]/=denominator
            if total_score<1:
                for GOterm,cscore in naive_GOterm_dict[Aspect]:
                    if GOterm in predict_dict:
                        predict_dict[GOterm]=total_score*predict_dict[GOterm
                            ]+(1-total_score)*cscore
                    else:
                        predict_dict[GOterm]=(1-total_score)*cscore
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
    annotation_dict=read_annotation()
    naive_GOterm_dict=read_naive_prob()
    target_list=read_fasta_as_list(infile)
    blast_dict=run_blast(infile)
    write_output(target_list,blast_dict,
        annotation_dict,naive_GOterm_dict,outfile)
