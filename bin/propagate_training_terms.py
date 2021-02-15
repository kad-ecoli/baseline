#!/usr/bin/env python
docstring='''propagate_training_term.py go-basic.obo uniprot_sprot_exp.txt
    propagate parent GO terms for CAFA training set "uniprot_sprot_exp.txt".
    output full sets of GO term to uniprot_sprot_exp.{F,P,C}
    output naive probability to naive.{F,P,C}
'''
import sys
import obo2csv # parsing GO hierachy
from math import log

def propagate_training_term(obo_dict,infile):
    allterm_dict={'F':dict(),'C':dict(),'P':dict()}
    fp=open(infile,'r')
    lines=fp.read().splitlines()
    fp.close()
    DB_Object_ID_list=[]
    for line in lines:
        DB_Object_ID,GO_ID,Aspect=line.split()
        if not DB_Object_ID in DB_Object_ID_list:
            DB_Object_ID_list.append(DB_Object_ID)
        if not GO_ID in obo_dict[Aspect]["Term"]:
            sys.stderr.write("ERROR! Cannot find GO Term %s\n"%GO_ID)
            continue
        if not DB_Object_ID in allterm_dict[Aspect]:
            allterm_dict[Aspect][DB_Object_ID]=[]
        if not GO_ID in allterm_dict[Aspect][DB_Object_ID]:
            allterm_dict[Aspect][DB_Object_ID].append(GO_ID)

        for parent_GO in obo_dict.is_a(Term_id=GO_ID, direct=False,
            name=True, number=False).split('\t'):
            GO_ID,name=parent_GO.split(" ! ")
            if not GO_ID in allterm_dict[Aspect][DB_Object_ID]:
                allterm_dict[Aspect][DB_Object_ID].append(GO_ID)

    for Aspect in "FPC":
        filename="uniprot_sprot_exp."+Aspect
        print("writing "+filename)
        txt=''
        for DB_Object_ID in DB_Object_ID_list:
            if not DB_Object_ID in allterm_dict[Aspect]:
                continue
            txt+="%s\t%s\n"%(DB_Object_ID,
                ','.join(sorted(allterm_dict[Aspect][DB_Object_ID])))
        fp=open(filename,'w')
        fp.write(txt)
        fp.close()

        filename="naive."+Aspect
        print("writing "+filename)
        GO_ID_list=[]
        for DB_Object_ID in allterm_dict[Aspect]:
            GO_ID_list+=allterm_dict[Aspect][DB_Object_ID]
        naive_list=[]
        ic_dict=dict()
        for GO_ID in sorted(set(GO_ID_list)):
            prob=1.*sum([GO_ID in allterm_dict[Aspect][DB_Object_ID] for \
                DB_Object_ID in allterm_dict[Aspect]])/len(allterm_dict[Aspect])
            name=obo_dict.short(GO_ID).split(' ! ')[1]
            ic_dict[GO_ID]=-log(prob,2)
            naive_list.append((GO_ID,prob,name))
        # we approximate the ic conditioned on all parent terms by the ic
        # conditioned on the most specific parent term
        ic_condition_on_parent=dict()
        for GO_ID,ic in ic_dict.items():
            ic_condition_on_parent[GO_ID]=ic
            max_parent_ic=0
            for parent_GO in obo_dict.is_a(Term_id=GO_ID, direct=False,
                name=True, number=False).split('\t'):
                parent_GO,name=parent_GO.split(" ! ")
                if parent_GO!=GO_ID and ic_dict[parent_GO]>max_parent_ic:
                    max_parent_ic=ic_dict[parent_GO]
            ic_condition_on_parent[GO_ID]-=max_parent_ic
        
        naive_list=[(GO_ID,prob,ic_dict[GO_ID],ic_condition_on_parent[GO_ID],name
            ) for GO_ID,prob,name in naive_list]
        naive_list.sort(key = lambda x: x[3])
        naive_list.sort(key = lambda x: x[1],reverse=True)

        txt=''
        for GO_ID,prob,ic,parent_ic,name in naive_list:
            txt+="%s\t%.6f\t%.6f\t%.6f\t%s"%(GO_ID,prob,ic,parent_ic,name)
        fp=open(filename,'w')
        fp.write(txt.replace("\t-0.000000","\t0.000000"))
        fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=3:
        sys.stderr.write(docstring)
        exit()

    obo_file=sys.argv[1]
    infile  =sys.argv[2]

    obo_dict=obo2csv.parse_obo_file(obo_file)
    propagate_training_term(obo_dict,infile)
