#!/usr/bin/env python
docstring='''predict_nw.py target.9606.fasta _9606_go.txt
    use nwalign to predict GO terms for input file target.9606.fasta,
    due to very long running time of nwalign, it will only be performed on
    templates found by a default blast run.
output:
    nwlocalID_1_9606_go.txt  - nident / length. blast baseline in CAFA
    nwalnscore_1_9606_go.txt - alnscore / self_alnscore
    nwglobalID_1_9606_go.txt - nident / qlen
    nwglobalID_2_9606_go.txt - nident / slen
    nwglobalID_3_9606_go.txt - nident / max(qlen,slen)
    nwrank_1_9606_go.txt     - 1 - rank / N
    nwfreq_1_9606_go.txt     - N(GOterm) / N
    nwmetago_1_9606_go.txt   - freq weighted by globalID, used in MetaGO
    nwnetgo_1_9606_go.txt    - freq weighted by alignment score, used in NetGO
'''
import os, sys
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

def read_fasta(infile):
    fasta_dict=dict()
    target_list=[]
    fp=open(infile,'r')
    for block in ('\n'+fp.read()).split('\n>')[1:]:
        lines=block.splitlines()
        target=lines[0].split()[0]
        sequence=''.join(lines[1:])
        target_list.append(target)
        fasta_dict[target]=sequence
    fp.close()
    return target_list,fasta_dict

def run_blast(infile):
    cmd="%s/blastp -db %s/uniprot_sprot_exp.fasta -outfmt '6 qacc sacc' -query %s"%(bindir,datdir,infile)
    p=Popen(cmd,shell=True,stdout=PIPE)
    stdout,stderr=p.communicate()
    blast_dict=dict()
    for line in stdout.splitlines():
        qacc,sacc=line.split('\t')
        if not qacc in blast_dict:
            blast_dict[qacc]=[]
        blast_dict[qacc].append(sacc)
    return blast_dict

def run_nw(target_list,fasta_dict,blast_dict,outfile):
    nw_dict=dict()
    for target in target_list:
        if not target in blast_dict:
            continue
        txt=">%s\n%s\n"%(target,fasta_dict[target])
        for suffix in [".target.fasta",".template.fasta"]:
            fp=open(outfile+suffix,'w')
            fp.write(txt)
            fp.close()

        cmd="%s/blastdbcmd -db %s/uniprot_sprot_exp.fasta -entry '%s' >> %s.template.fasta"%(
            bindir,datdir,','.join(blast_dict[target]),outfile)
        os.system(cmd)
        
        cmd="%s/NWalign %s.target.fasta %s.template.fasta | sort -k8nr"%(
            bindir,outfile,outfile)
        p=Popen(cmd,shell=True,stdout=PIPE)
        stdout,stderr=p.communicate()

        self_alnscore=1
        nw_dict[target]=[]
        for line in stdout.splitlines():
            if line.startswith('#'):
                continue
            target,template,ID1,ID2,IDali,L1,L2,Lali,score=line.split('\t')
            if target==template:
                self_alnscore=float(score)
                continue
            IDali=float(IDali)
            ID1=float(ID1)
            ID2=float(ID2)
            score=max(0,float(score))
            nw_dict[target].append([template,
                IDali,                    # localID
                score,                    # alnscore / self_alnscore
                ID1,                      # globalID1
                ID2,                      # globalID2
                min((ID1,ID2)),           # globalID3
                -1.*len(nw_dict[target]), # ranking
                1.,                       # freq
                ID1,                      # MetaGO
                score,                    # NetGO
            ])
        for t in range(len(nw_dict[target])):
            nw_dict[target][t][2]=min(1,nw_dict[target][t][2]/self_alnscore)

    cmd="rm %s.target.fasta %s.template.fasta"%(outfile,outfile)
    os.system(cmd)
    return nw_dict

def write_output(nw_dict,annotation_dict,suffix):
    method_list=[
        "nwlocalID_1",  # 0 - nident / length. blast baseline in CAFA
        "nwalnscore_1", # 1 - alnscore / self_alnscore
        "nwglobalID_1", # 2 - nident / qlen
        "nwglobalID_2", # 3 - nident / slen
        "nwglobalID_3", # 4 - nident / max(qlen,slen)
        "nwrank_1",     # 5 - 1 - rank / N
        "nwfreq_1",     # 6 - N(GOterm) / N
        "nwmetago_1",   # 7 - freq weighted by globalID, used in MetaGO
        "nwnetgo_1",    # 8 - freq weighted by bitscore, used in NetGO
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
                            ]/len(nw_dict[target])
                elif m>=6:
                    if denominator<=0:
                        continue
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
    target_list,fasta_dict=read_fasta(infile)
    blast_dict=run_blast(infile)
    nw_dict=run_nw(target_list,fasta_dict,blast_dict,suffix)
    write_output(nw_dict,annotation_dict,suffix)
