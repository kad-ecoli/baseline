#!/usr/bin/env python
docstring='''predict_blast.py target.9606.fasta _9606_go.txt
    use blast baseline to predict GO terms for input file target.9606.fasta,
output:
    blastlocalID_1_9606_go.txt  - nident / length. blast baseline in CAFA
    blastglobalID_1_9606_go.txt - nident / qlen
    blastglobalID_2_9606_go.txt - nident / slen
    blastglobalID_3_9606_go.txt - nident / max(qlen,slen)
    blastevalue_1_9606_go.txt   - 1-sigmoid(evalue)
    blastrank_1_9606_go.txt     - 1 - rank / N
    blastfreq_1_9606_go.txt     - N(GOterm) / N
    blastmetago_1_9606_go.txt   - freq weighted by globalID, used in MetaGO
    blastnetgo_1_9606_go.txt    - freq weighted by bitscore, used in NetGO
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

def run_blast(infile):
    cmd="%s/blastp -db %s/uniprot_sprot_exp.fasta -outfmt '6 qacc qlen sacc slen evalue bitscore length nident' -query %s"%(bindir,datdir,infile)
    p=Popen(cmd,shell=True,stdout=PIPE)
    stdout,stderr=p.communicate()
    blast_dict=dict()
    for line in stdout.splitlines():
        qacc,qlen,sacc,slen,evalue,bitscore,length,nident=line.split('\t')
        qlen=float(qlen)
        slen=float(slen)
        evalue=float(evalue)
        bitscore=float(bitscore)
        length=float(length)
        nident=float(nident)
        if not qacc in blast_dict:
            blast_dict[qacc]=[]
        blast_dict[qacc].append([sacc,
            nident/length,             # localID
            nident/qlen,               # globalID1
            nident/slen,               # globalID2
            nident/max(qlen,slen),     # globalID3
            1-1/(1+exp(-evalue)),      # evalue
            -1.*len(blast_dict[qacc]), # ranking
            1.,                        # freq
            nident/qlen,               # MetaGO
            bitscore,                  # NetGO
            ])
    return blast_dict

def write_output(blast_dict,annotation_dict,suffix):
    method_list=[
        "blastlocalID_1",  # 0 - nident / length. blast baseline in CAFA
        "blastglobalID_1", # 1 - nident / qlen
        "blastglobalID_2", # 2 - nident / slen
        "blastglobalID_3", # 3 - nident / max(qlen,slen)
        "blastevalue_1",   # 4 - 1 - sigmoid(evalue)
        "blastrank_1",     # 5 - 1 - rank / N
        "blastfreq_1",     # 5 - N(GOterm) / N
        "blastmetago_1",   # 6 - freq weighted by globalID, used in MetaGO
        "blastnetgo_1",    # 7 - freq weighted by bitscore, used in NetGO
    ]

    for m,method in enumerate(method_list):
        outfile=method+suffix
        print("writing "+outfile)
        author,model=method.split('_')
        txt ="AUTHOR %s\n"%author
        txt+="MODEL %s\n"%model
        txt+="KEYWORDS sequence alignment.\n"
        for target in sorted(blast_dict.keys()):
            for Aspect in "FPC":
                predict_dict=dict()
                denominator=0
                for line in blast_dict[target]:
                    sacc=line[0]
                    score=line[m+1]
                    if not sacc in annotation_dict[Aspect]:
                        continue
                    GOterm_list=annotation_dict[Aspect][sacc]
                    if m<=5:
                        for GOterm in GOterm_list:
                            if not GOterm in predict_dict or \
                                score>predict_dict[GOterm]:
                                predict_dict[GOterm]=score
                    elif m>=6:
                        denominator+=score
                        for GOterm in GOterm_list:
                            if not GOterm in predict_dict:
                                predict_dict[GOterm]=0
                            predict_dict[GOterm]+=score
                if m==5:
                    for GOterm in predict_dict:
                        predict_dict[GOterm]=1+predict_dict[GOterm
                            ]/len(blast_dict[target])
                elif m>=6:
                    for GOterm in predict_dict:
                        predict_dict[GOterm]/=denominator
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
    annotation_dict=read_annotation()
    blast_dict=run_blast(infile)
    write_output(blast_dict,annotation_dict,suffix)
