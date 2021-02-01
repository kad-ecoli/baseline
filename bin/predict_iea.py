#!/usr/bin/env python
docstring='''predict_iea.py target.9606.fasta iea_1_9606_go.txt
    use IEA annotation to predict GO terms for input file target.9606.fasta,
    output predictions to naive_1_9606_go.txt
'''
import sys
from os.path import dirname, basename, abspath
bindir=dirname(abspath(__file__))
rootdir=dirname(bindir)
datdir=rootdir+"/data"

def read_precision_by_evidence():
    evidence2cscore_dict=dict()
    fp=open(datdir+"/precision_by_evidence.txt",'r')
    for line in fp.read().splitlines():
        if line.startswith('#'):
            continue
        evidence,cscore=line.split()
        evidence2cscore_dict[evidence]=float(cscore)
    fp.close()
    return evidence2cscore_dict

def read_target_map():
    entry_dict=dict()
    fp=open(rootdir+"/input/target.map",'r')
    for line in fp.read().splitlines():
        entry,accession=line.split()
        entry_dict[entry]=accession
    return entry_dict

def read_fasta_as_list(infile,entry_dict):
    target_list=[]
    target_dict=dict()
    fp=open(infile,'r')
    for line in fp.read().splitlines():
        if line.startswith('>'):
            target,entry=line[1:].split()
            if not entry in entry_dict:
                continue
            target_list.append(target)
            accession=entry_dict[entry]
            target_dict[target]=accession
    fp.close()
    return target_list, target_dict

def read_goa_isa(target_dict,evidence2cscore_dict):
    accession_set=set(target_dict.values())
    annotation_dict=dict()
    fp=open(datdir+"/goa_uniprot_all.is_a",'r')
    for line in fp.read().splitlines():
        accession,Aspect,GOterm,evidence_list=line.split('\t')
        if not accession in accession_set:
            continue
        cscore=0
        for evidence in evidence_list.split(','):
            if evidence2cscore_dict[evidence]>cscore:
                cscore=evidence2cscore_dict[evidence]
        if not accession in annotation_dict:
            annotation_dict[accession]=[]
        annotation_dict[accession].append((cscore,GOterm))
    fp.close()
    return annotation_dict

def write_output(target_list,target_dict,annotation_dict,outfile):
    txt ="AUTHOR iea\n"
    txt+="MODEL 1\n"
    txt+="KEYWORDS de novo prediction.\n"
    for target in target_list:
        accession=target_dict[target]
        if not accession in annotation_dict:
            continue
        for cscore,GOterm in annotation_dict[accession]:
            txt+="%s\t%s\t%.2f\n"%(target,GOterm,cscore)
    txt+="END\n"
    fp=open(outfile,'w')
    fp.write(txt)
    return

if __name__=="__main__":
    if len(sys.argv)!=3:
        sys.stderr.write(docstring)
        exit()

    infile                    =sys.argv[1]
    outfile                   =sys.argv[2]

    evidence2cscore_dict=read_precision_by_evidence()
    entry_dict=read_target_map()
    target_list, target_dict=read_fasta_as_list(infile,entry_dict)
    annotation_dict=read_goa_isa(target_dict,evidence2cscore_dict)
    write_output(target_list,target_dict,annotation_dict,outfile)
