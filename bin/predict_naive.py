#!/usr/bin/env python
docstring='''predict_naive.py target.9606.fasta _9606_go.txt 9606
    use naive baseline to predict GO terms for input file target.9606.fasta,
Output:
    naive_1_9606_go.txt - naive trained on the whole database
    naive_2_9606_go.txt - naive trained on the same species proteins
                          from the database
'''
import sys
from os.path import dirname, basename, abspath, isfile
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

def train_species_naive_prob(species):
    species_naive_GOterm_list=[]
    template_list=[]
    fp=open(datdir+"/uniprot_sprot_exp.species",'r')
    for line in fp.read().splitlines():
        items=line.split()
        if items[1]==species:
            template_list.append(items[0])
    fp.close()
    sys.stderr.write("training on %d templates\n"%len(template_list))
    template_set=set(template_list)

    for Aspect in "FPC":
        annotation_dict=dict()
        GOterm_list=[]
        fp=open(datdir+"/uniprot_sprot_exp."+Aspect,'r')
        for line in fp.read().splitlines():
            template,GOterms=line.split()
            if template in template_set:
                annotation_dict[template]=GOterms.split(',')
                GOterm_list+=annotation_dict[template]
                annotation_dict[template]=set(annotation_dict[template])
        fp.close()

        if len(annotation_dict)==0:
            filename=datdir+"/naive."+Aspect
            sys.stderr.write("0 template for %s aspect. fallback to %s\n"%(
                Aspect,filename))
            fp=open(filename,'r')
            for line in fp.read().splitlines():
                GOterm,cscore=line.split('\t')[:2]
                if "%.2f"%float(cscore)=="0.00":
                    break
                species_naive_GOterm_list.append((float(cscore),GOterm))
            fp.close()
        else:
            cscore_list=[]
            for GOterm in list(set(GOterm_list)):
                cscore=(1.+sum([GOterm in annotation_dict[template] for \
                    template in annotation_dict]))/(1.+len(annotation_dict))
                if "%.2f"%float(cscore)=="0.00":
                    continue
                cscore_list.append((cscore,GOterm))
            cscore_list.sort(reverse=True)
            species_naive_GOterm_list+=cscore_list
    return sorted(species_naive_GOterm_list,reverse=True)

def read_fasta_as_list(infile):
    target_list=[]
    fp=open(infile,'r')
    for line in fp.read().splitlines():
        if line.startswith('>'):
            target_list.append(line[1:].split()[0])
    fp.close()
    return target_list

def write_output(target_list,naive_GOterm_list,species_naive_GOterm_list,suffix):
    txt ="AUTHOR naive\n"
    txt+="MODEL 1\n"
    txt+="KEYWORDS de novo prediction.\n"
    for target in target_list:
        for cscore,GOterm in naive_GOterm_list:
            txt+="%s\t%s\t%.2f\n"%(target,GOterm,cscore)
    txt+="END\n"
    fp=open("naive_1"+suffix,'w')
    fp.write(txt)
    fp.close()

    if len(species_naive_GOterm_list)==0:
        return
    txt ="AUTHOR naive\n"
    txt+="MODEL 2\n"
    txt+="KEYWORDS de novo prediction.\n"
    for target in target_list:
        for cscore,GOterm in species_naive_GOterm_list:
            txt+="%s\t%s\t%.2f\n"%(target,GOterm,cscore)
    txt+="END\n"
    fp=open("naive_2"+suffix,'w')
    fp.write(txt)
    fp.close()
    return

if __name__=="__main__":
    if not len(sys.argv) in [3,4]:
        sys.stderr.write(docstring)
        exit()

    infile =sys.argv[1]
    suffix =sys.argv[2]
    species=''
    species_naive_GOterm_list=[]
    if len(sys.argv)>3:
        species=sys.argv[3]
        species_naive_GOterm_list=train_species_naive_prob(species)
    naive_GOterm_list=read_naive_prob()
    target_list=read_fasta_as_list(infile)
    write_output(target_list,naive_GOterm_list,species_naive_GOterm_list,suffix)
