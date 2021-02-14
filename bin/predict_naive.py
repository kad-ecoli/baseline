#!/usr/bin/env python
docstring='''predict_naive.py target.9606.fasta naive_1_9606_go.txt
    use naive baseline to predict GO terms for input file target.9606.fasta,
    output predictions to naive_1_9606_go.txt
'''
import sys
from os.path import dirname, basename, abspath
bindir=dirname(abspath(__file__))
rootdir=dirname(bindir)
datdir=rootdir+"/data"

def read_naive_prob():
    naive_GOterm_list=[]
    for Aspect in "FPC":
        filename=datdir+"/naive."+Aspect
        fp=open(filename,'r')
        for line in fp.read().splitlines():
            GOterm,cscore=line.split('\t')[:2]
            if "%.2f"%float(cscore)=="0.00":
                break
            naive_GOterm_list.append((float(cscore),GOterm))
        fp.close()
    return sorted(naive_GOterm_list,reverse=True)

def read_fasta_as_list(infile):
    target_list=[]
    fp=open(infile,'r')
    for line in fp.read().splitlines():
        if line.startswith('>'):
            target_list.append(line[1:].split()[0])
    fp.close()
    return target_list

def write_output(target_list,naive_GOterm_list,outfile):
    txt ="AUTHOR naive\n"
    txt+="MODEL 1\n"
    txt+="KEYWORDS de novo prediction.\n"
    for target in target_list:
        for cscore,GOterm in naive_GOterm_list:
            txt+="%s\t%s\t%.2f\n"%(target,GOterm,cscore)
    txt+="END\n"
    fp=open(outfile,'w')
    fp.write(txt)
    return

if __name__=="__main__":
    if len(sys.argv)!=3:
        sys.stderr.write(docstring)
        exit()

    infile =sys.argv[1]
    outfile=sys.argv[2]
    naive_GOterm_list=read_naive_prob()
    target_list=read_fasta_as_list(infile)
    write_output(target_list,naive_GOterm_list,outfile)
