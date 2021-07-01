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
root_terms={ "GO:0005575", "GO:0003674", "GO:0008150"}

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
            GOterms=set([GOterm for GOterm in GOterms.split(','
                           ) if not GOterm in root_terms])
            if target in NK_list:
                label_dict['NK'][Aspect][target]=GOterms
            elif target in LK_list:
                label_dict['LK'][Aspect][target]=GOterms
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
        if GOterm in root_terms:
            continue
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

def sum_ic(GOterm_list,ic_dict):
    return sum([ic_dict[GOterm] for GOterm in GOterm_list if GOterm in ic_dict])
    

def assess_result(label_dict,ic_dict,predict_dict,outfile):
    txt="Aspect\tType\tMode\tFmax\tCutoff\tSmin\tCutoff\twFmax\tCutoff\tCoverage\n"
    aspect_dict=dict(F="mfo",P="bpo",C="cco")
    for Aspect in "FPC":
        ontology=aspect_dict[Aspect]
        for target_type in ["NK","LK"]:
            Fmax=0
            Cutoff1=0
            Smin=0
            Cutoff2=1
            wFmax=0
            Cutoff3=0
            total_label=len(label_dict[target_type][Aspect])
            Coverage=1.*len(predict_dict[target_type][Aspect])/total_label
            cscore_list=[]
            label_ic_dict=dict()
            for target in label_dict[target_type][Aspect]:
                label_ic_dict[target]=sum_ic(
                    label_dict[target_type][Aspect][target],ic_dict)
                Smin+=label_ic_dict[target]
            Smin/=total_label
            for target,GOterms in predict_dict[target_type][Aspect].items():
                cscore_list+=[cscore for GOterm,cscore in GOterms]
            cscore_list=sorted(set(cscore_list))
            for cutoff in cscore_list:
                total_precision=0
                total_recall=0
                total_predict=0
                total_wprecision=0
                total_wrecall=0
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
                    wpredict=sum_ic(predict_list,ic_dict)
                    tp=label_set.intersection(predict_list)
                    wtp=sum_ic(tp,ic_dict)
                    precision=1.*len(tp)/len(predict_list)
                    wprecision=0 if wpredict==0 else wtp/wpredict
                    recall=1.*len(tp)/len(label_set)
                    wrecall=0 if label_ic_dict[target]==0 else wtp/label_ic_dict[target]
                    total_precision+=precision
                    total_recall+=recall
                    total_wprecision+=wprecision
                    total_wrecall+=wrecall
                    for GOterm in predict_list:
                        if GOterm in ic_dict and not GOterm in label_set:
                            total_mi+=ic_dict[GOterm]
                    for GOterm in label_set:
                        if GOterm in ic_dict and not GOterm in predict_list:
                            total_ru+=ic_dict[GOterm]
                total_label=len(label_dict[target_type][Aspect])
                precision=0 if total_predict==0 else total_precision/total_predict
                wprecision=0 if total_predict==0 else total_wprecision/total_predict
                recall   =total_recall/total_label
                wrecall  =total_wrecall/total_label
                mi       =total_mi/total_label
                ru       =total_ru/total_label

                F=0 if precision*recall==0 else 2/(1/precision+1/recall)
                wF=0 if wprecision*wrecall==0 else 2/(1/wprecision+1/wrecall)
                S=sqrt(ru*ru+mi*mi)
                if F>=Fmax:
                    Fmax=F
                    Cutoff1=cutoff
                if S<=Smin:
                    Smin=S
                    Cutoff2=cutoff
                if wF>=wFmax:
                    wFmax=wF
                    Cutoff3=cutoff
            txt+="%s\t%s\tfull\t%.4f\t%.2f\t%.4f\t%.2f\t%.4f\t%.2f\t%.4f\n"%(
                ontology,target_type.upper(),
                Fmax,Cutoff1,Smin,Cutoff2,wFmax,Cutoff3,Coverage)
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
