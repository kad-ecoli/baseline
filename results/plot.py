#!/usr/bin/env python
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import os
from os.path import dirname, basename, abspath
resultdir=dirname(abspath(__file__))
rootdir=dirname(resultdir)

method_list=[
    "blastevalue_1",
    "blastlocalID_1",
    "naive_1",
    "naive_2",
    "blastbitscore_1",
    "blastbitscore_2",
    "blastbitscore_3",
    "blastglobalID_1",
    "blastglobalID_2",
    "blastglobalID_3",
    "blastrank_1",
    "iea_1",
    "blastfreq_1",
    "blastmetago_1",
    "blastnetgo_1",
    ]
method_list=[method for method in method_list if os.path.isfile(
    "%s/%s_all_results.txt"%(resultdir,method))]
fontsize=9

for metric in ["Fmax","Smin","wFmax"]:
    plt.figure(figsize=(7.87,7.87))
    for a,Aspect in enumerate(["mf","bp","cc"]):
        for k,Knowledge in enumerate(["NK","LK"]):
            ax=plt.subplot(3,2,2*a+k+1)

            score_list=[]
            for m,method in enumerate(method_list):
                infile="%s/%s_all_results.txt"%(resultdir,method)
                score=0
                if os.path.isfile(infile):
                    fp=open(infile,'r')
                    for line in fp.read().splitlines():
                        if line.startswith("%so\t%s\tfull"%(Aspect,Knowledge)):
                            if metric=="Fmax":
                                score=float(line.split()[3])
                            elif metric=="Smin":
                                score=float(line.split()[5])
                            elif metric=="wFmax":
                                score=float(line.split()[7])
                    fp.close()
                score_list.append(score)
                plt.text(m,score+0.01,("%.3f"%score).lstrip('0')[:4],
                    rotation=90,fontsize=fontsize,va="bottom",ha="center")
            plt.bar(range(len(score_list)),score_list,
                align="center",color="lightgrey")
            plt.ylabel("%s %s (%s targets)"%(Aspect.upper(),metric,Knowledge
                 ), labelpad=0,fontsize=fontsize)
            plt.yticks(fontsize=fontsize)
            plt.xticks(range(len(score_list)),
                [method.replace('_','') for method in method_list],
                rotation=90,fontsize=fontsize)
            ymax=0.8
            if metric=="Smin":
                ymax=20
            plt.axis([-0.5,len(method_list)-0.5,0,ymax])
            ax.tick_params('x',length=0)
            ax.tick_params('y',length=3)
            ax.tick_params('both',pad=0)
    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    plt.savefig(metric+"_full.png",dpi=300)
