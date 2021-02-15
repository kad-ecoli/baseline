#!/usr/bin/env python
docstring='''propagate_groundtruth_term.py go-basic.obo F leafonly_MFO.txt leafonly_MFO.is_a
    propagate parent GO terms for CAFA test set leafonly_MFO.txt
    output full sets of GO term to leafonly_MFO.is_a
'''
import sys
import obo2csv # parsing GO hierachy
from math import log

def propagate_groundtruth_term(obo_dict,Aspect,infile,outfile):
    allterm_dict=dict()
    fp=open(infile,'r')
    lines=fp.read().splitlines()
    fp.close()
    DB_Object_ID_list=[]
    missing_GO_list=[]
    for line in lines:
        DB_Object_ID,GO_ID=line.split()
        if not DB_Object_ID in DB_Object_ID_list:
            DB_Object_ID_list.append(DB_Object_ID)
        if not GO_ID in obo_dict[Aspect]["Term"]:
            if not GO_ID in missing_GO_list:
                sys.stderr.write("ERROR! Cannot find GO Term %s\n"%GO_ID)
                missing_GO_list.append(GO_ID)
            continue
        if not DB_Object_ID in allterm_dict:
            allterm_dict[DB_Object_ID]=[]
        if not GO_ID in allterm_dict[DB_Object_ID]:
            allterm_dict[DB_Object_ID].append(GO_ID)

        for parent_GO in obo_dict.is_a(Term_id=GO_ID, direct=False,
            name=True, number=False).split('\t'):
            GO_ID,name=parent_GO.split(" ! ")
            if not GO_ID in allterm_dict[DB_Object_ID]:
                allterm_dict[DB_Object_ID].append(GO_ID)

    print("writing "+outfile)
    txt=''
    for DB_Object_ID in DB_Object_ID_list:
        if not DB_Object_ID in allterm_dict:
            continue
        txt+="%s\t%s\n"%(DB_Object_ID,
            ','.join(sorted(allterm_dict[DB_Object_ID])))
    fp=open(outfile,'w')
    fp.write(txt)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=5:
        sys.stderr.write(docstring)
        exit()

    obo_file=sys.argv[1]
    Aspect  =sys.argv[2]
    infile  =sys.argv[3]
    outfile =sys.argv[4]

    obo_dict=obo2csv.parse_obo_file(obo_file)
    propagate_groundtruth_term(obo_dict,Aspect,infile,outfile)
