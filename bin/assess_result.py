#!/usr/bin/env python
docstring='''assess_result.py blastlocalID_1_all.txt blastlocalID_1_all_results.txt
    calculate Fmax, coverage and Smin by full mode for prediction 
    blastlocalID_1_all.txt. output result to blastlocalID_1_all_results.txt
'''
import sys
import obo2csv # parsing GO hierachy
from os.path import dirname, basename, abspath
from math import sqrt
bindir=dirname(abspath(__file__))
rootdir=dirname(bindir)
datdir=rootdir+"/data"
labeldir=rootdir+"/groundtruth"

def read_list_file(list_file):
    fp=open(list_file,'r')
    target_list=fp.read().splitlines()
    fp.close()
    return target_list

def read_label():
    label_dict=dict(NK=dict(F=dict(),P=dict(),C=dict()),
                    LK=dict(F=dict(),P=dict(),C=dict()))
    NK_list =read_list_file(labeldir+"/all_type1.txt")
    LK_list =read_list_file(labeldir+"/all_type2.txt")
    aspect_dict=dict(F="mfo",P="bpo",C="cco")
    for Aspect in "FPC":
        filename=labeldir+"/leafonly_"+aspect_dict[Aspect].upper()+".is_a"
        fp=open(filename,'r')
        lines=fp.read().splitlines()
        fp.close()
        for line in lines:
            target,GOterms=line.split()
            if target in NK_list:
                label_dict['NK'][Aspect][target]=set(GOterms.split(','))
            elif target in LK_list:
                label_dict['LK'][Aspect][target]=set(GOterms.split(','))
            else:
                print("WARNING! Unclassified labeled target "+target)
    for target in NK_list:
        num_aspect=sum([target in label_dict["NK"][Aspect] for Aspect in "FPC"])
        if num_aspect==0:
            print("WARNING! 0 label for NK target "+target)
    for target in LK_list:
        num_aspect=sum([target in label_dict["LK"][Aspect] for Aspect in "FPC"])
        if num_aspect==0:
            print("WARNING! 0 label for LK target "+target)
    return label_dict

def read_information_content():
    ic_dict=dict()
    for Aspect in "FPC":
        filename=datdir+"/naive."+Aspect
        fp=open(filename,'r')
        for line in fp.read().splitlines():
            items=line.split('\t')
            GOterm=items[0]
            ic=float(items[3])
            ic_dict[GOterm]=ic
        fp.close()
    return ic_dict

def read_prediction(obo_dict,infile,label_dict):
    predict_dict=dict(NK=dict(F=dict(),P=dict(),C=dict()),
                      LK=dict(F=dict(),P=dict(),C=dict()))
    fp=open(infile,'r')
    lines=fp.read().splitlines()
    fp.close()
    for line in lines:
        items=line.split('\t')
        if len(items)!=3:
            continue
        target,GOterm,cscore=items
        Aspect=''
        for a in "FPC":
            if GOterm in obo_dict[a]["Term"]:
                Aspect=a
                break
        if Aspect=='':
            print("ERROR! Unknow GO term "+GOterm)
            continue
        target_type=''
        for t in ["NK","LK"]:
            if target in label_dict[t][Aspect]:
                target_type=t
        if target_type=='':
            continue
        if not target in predict_dict[target_type][Aspect]:
            predict_dict[target_type][Aspect][target]=[]
        predict_dict[target_type][Aspect][target].append((GOterm,float(cscore)))
    return predict_dict

def assess_result(label_dict,ic_dict,predict_dict,outfile):
    txt="Aspect\tType\tMode\tFmax\tCutoff\tSmin\tCutoff\tCoverage\n"
    aspect_dict=dict(F="mfo",P="bpo",C="cco")
    for Aspect in "FPC":
        ontology=aspect_dict[Aspect]
        for target_type in ["NK","LK"]:
            Fmax=0
            Cutoff1=0
            Smin=0
            Cutoff2=1
            total_label=len(label_dict[target_type][Aspect])
            Coverage=1.*len(predict_dict[target_type][Aspect])/total_label
            cscore_list=[]
            label_ic_dict=dict()
            for target in label_dict[target_type][Aspect]:
                label_ic_dict[target]=0
                for GOterm in label_dict[target_type][Aspect][target]:
                    if GOterm in ic_dict:
                        label_ic_dict[target]+=ic_dict[GOterm]
                Smin+=label_ic_dict[target]
            Smin/=total_label
            for target,GOterms in predict_dict[target_type][Aspect].items():
                cscore_list+=[cscore for GOterm,cscore in GOterms]
            cscore_list=sorted(set(cscore_list))
            for cutoff in cscore_list:
                total_precision=0
                total_recall=0
                total_predict=0
                total_ru=0
                total_mi=0
                for target in label_dict[target_type][Aspect]:
                    if not target in predict_dict[target_type][Aspect]:
                        total_ru+=label_ic_dict[target]
                        continue
                    predict_list=[GOterm for GOterm,cscore in predict_dict[
                        target_type][Aspect][target] if cscore>=cutoff]
                    if len(predict_list)==0:
                        total_ru+=label_ic_dict[target]
                        continue
                    label_set=label_dict[target_type][Aspect][target]
                    total_predict+=1
                    tp=label_set.intersection(predict_list)
                    precision=1.*len(tp)/len(predict_list)
                    recall=1.*len(tp)/len(label_set)
                    total_precision+=precision
                    total_recall+=recall
                    for GOterm in predict_list:
                        if GOterm in ic_dict and not GOterm in label_set:
                            total_mi+=ic_dict[GOterm]
                    for GOterm in label_set:
                        if GOterm in ic_dict and not GOterm in predict_list:
                            total_ru+=ic_dict[GOterm]
                total_label=len(label_dict[target_type][Aspect])
                precision=total_precision/total_predict
                recall   =total_recall/total_label
                mi       =total_mi/total_label
                ru       =total_ru/total_label

                F=2/(1/precision+1/recall)
                S=sqrt(ru*ru+mi*mi)
                if F>=Fmax:
                    Fmax=F
                    Cutoff1=cutoff
                if S<=Smin:
                    Smin=S
                    Cutoff2=cutoff
            txt+="%s\t%s\tfull\t%.4f\t%.2f\t%.4f\t%.2f\t%.4f\n"%(
                ontology,target_type.upper(),Fmax,Cutoff1,Smin,Cutoff2,Coverage)
    fp=open(outfile,'w')
    fp.write(txt)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=3:
        sys.stderr.write(docstring)
        exit()

    infile  =sys.argv[1]
    outfile =sys.argv[2]

    obo_dict=obo2csv.parse_obo_file(datdir+"/go-basic.obo")
    label_dict=read_label()
    ic_dict=read_information_content()
    predict_dict=read_prediction(obo_dict,infile,label_dict)
    assess_result(label_dict,ic_dict,predict_dict,outfile)
