#!/usr/bin/env python
docstring='''predict_nw.py target.9606.fasta _9606_go.txt
    use nwalign to predict GO terms for input file target.9606.fasta,
    for consistency with blast, only top 500 hits with largest alignment
    scores are considered.
output:
    nwlocalID_1_9606_go.txt  - nident / length. blast baseline in CAFA
    nwglobalID_1_9606_go.txt - nident / qlen
    nwglobalID_2_9606_go.txt - nident / slen
    nwglobalID_3_9606_go.txt - nident / max(qlen,slen)
    nwrank_1_9606_go.txt     - 1 - rank / N
    nwfreq_1_9606_go.txt     - N(GOterm) / N
    nwmetago_1_9606_go.txt   - freq weighted by globalID, used in MetaGO
    nwnetgo_1_9606_go.txt    - freq weighted by alignment score, used in NetGO
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

def run_nw(infile):
    cmd="%s/NWalign %s %s/uniprot_sprot_exp.fasta|sort -k9nr"%(bindir,infile,datdir)
    p=Popen(cmd,shell=True,stdout=PIPE)
    stdout,stderr=p.communicate()
    nw_dict=dict()
    for line in stdout.splitlines():
        if line.startswith('#'):
            continue
        target,template,ID1,ID2,IDali,L1,L2,Lali,score=line.split('\t')
        if not target in nw_dict:
            nw_dict[target]=[]
        elif len(nw_dict[target])>=500:
            continue
        IDali=float(IDali)
        ID1=float(ID1)
        ID2=float(ID2)
        score=float(score)
        nw_dict[target].append([template,
            IDali,                    # localID
            ID1,                      # globalID1
            ID2,                      # globalID2
            min((ID1,ID2)),           # globalID3
            -1.*len(nw_dict[target]), # ranking
            1.,                       # freq
            ID1,                      # MetaGO
            score,                    # NetGO
            ])
    return nw_dict

def write_output(nw_dict,annotation_dict,suffix):
    method_list=[
        "nwlocalID_1",  # 0 - nident / length. blast baseline in CAFA
        "nwglobalID_1", # 1 - nident / qlen
        "nwglobalID_2", # 2 - nident / slen
        "nwglobalID_3", # 3 - nident / max(qlen,slen)
        "nwrank_1",     # 4 - 1 - rank / N
        "nwfreq_1",     # 5 - N(GOterm) / N
        "nwmetago_1",   # 6 - freq weighted by globalID, used in MetaGO
        "nwnetgo_1",    # 7 - freq weighted by bitscore, used in NetGO
    ]

    for m,method in enumerate(method_list):
        outfile=method+suffix
        print("writing "+outfile)
        author,model=method.split('_')
        txt ="AUTHOR %s\n"%author
        txt+="MODEL %s\n"%model
        txt+="KEYWORDS sequence alignment.\n"
        for target in sorted(nw_dict.keys()):
            for Aspect in "FPC":
                predict_dict=dict()
                denominator=0
                for line in nw_dict[target]:
                    template=line[0]
                    score=line[m+1]
                    if not template in annotation_dict[Aspect]:
                        continue
                    GOterm_list=annotation_dict[Aspect][template]
                    if m<=4:
                        for GOterm in GOterm_list:
                            if not GOterm in predict_dict or \
                                score>predict_dict[GOterm]:
                                predict_dict[GOterm]=score
                    elif m>=5:
                        denominator+=score
                        for GOterm in GOterm_list:
                            if not GOterm in predict_dict:
                                predict_dict[GOterm]=0
                            predict_dict[GOterm]+=score
                if m==4:
                    for GOterm in predict_dict:
                        predict_dict[GOterm]=1+predict_dict[GOterm
                            ]/len(nw_dict[target])
                elif m>=5:
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
    nw_dict=run_nw(infile)
    write_output(nw_dict,annotation_dict,suffix)
