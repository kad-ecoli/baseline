#!/usr/bin/env python
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import os
from os.path import dirname, basename, abspath
resultdir=dirname(abspath(__file__))
rootdir=dirname(resultdir)

method_list=[
    "localID_1",
    "alnscore_1",
    "alnscore_2",
    "alnscore_3",
    "globalID_1",
    "globalID_2",
    "globalID_3",
    "rank_1",
    "freq_1",
    "metago_1",
    "netgo_1",
    ]
method_list=[method for method in method_list if os.path.isfile(
    "%s/blast%s_all_results.txt"%(resultdir,method.replace("alnscore","bitscore")))
    or os.path.isfile(
    "%s/nw%s_all_results.txt"%(resultdir,method))
    ]
fontsize=9
color_list=["darkgrey","lightgrey"]

for metric in ["Fmax","Smin","wFmax"]:
    plt.figure(figsize=(7.87,7.87))
    for a,Aspect in enumerate(["mf","bp","cc"]):
        for k,Knowledge in enumerate(["NK","LK"]):
            ax=plt.subplot(3,2,2*a+k+1)

            for t,tool in enumerate(["nw","blast"]):
                score_list=[]
                for m,method in enumerate(method_list):
                    if t==1 and method.startswith("alnscore"):
                        method=method.replace("alnscore","bitscore")
                    infile="%s/%s%s_all_results.txt"%(
                        resultdir,tool,method)
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
                        rotation=90, fontsize=fontsize,va="bottom",
                        ha="right" if t==0 else "left")
                plt.bar(range(len(score_list)),score_list, width=0.4*(2*t-1),
                    align="edge", color=color_list[t], label=tool)
            plt.ylabel("%s %s (%s targets)"%(Aspect.upper(),metric,Knowledge),
                labelpad=0,fontsize=fontsize)
            plt.yticks(fontsize=fontsize)
            xticks=[]
            for method in method_list:
                xticks.append(method.replace('_',''))
                for t in "123":
                    if xticks[-1]=="alnscore"+t:
                        xticks[-1]+="\nbitscore"+t
            plt.xticks(range(len(score_list)),xticks,rotation=90,fontsize=fontsize)
            ymax=0.8
            if metric=="Smin":
                ymax=20
            plt.axis([-0.5,len(method_list)-0.5,0,ymax])
            ax.tick_params('x',length=0)
            ax.tick_params('y',length=3)
            ax.tick_params('both',pad=0)
            if a==0 and k==1:
                plt.legend(fontsize=fontsize, ncol=2, loc="upper left",
                    borderpad=0.1,handlelength=1,handletextpad=0.1,
                    borderaxespad=0.2,columnspacing=1)
    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    plt.savefig(metric+"_full_nw.png",dpi=300)
