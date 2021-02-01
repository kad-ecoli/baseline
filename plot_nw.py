#!/usr/bin/env python
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import os
from os.path import dirname, basename, abspath
rootdir=dirname(abspath(__file__))

method_list=[
    "localID_1",
    "alnscore_1",
    "globalID_1",
    "globalID_2",
    "globalID_3",
    "rank_1",
    "freq_1",
    "metago_1",
    "netgo_1",
    ]
method_list=[method for method in method_list if os.path.isfile(
    "%s/CAFA_assessment_tool/results/blast%s_all_results.txt"%(rootdir,method))
    or os.path.isfile(
    "%s/CAFA_assessment_tool/results/nw%s_all_results.txt"%(rootdir,method))
    ]
fontsize=9
color_list=["darkgrey","lightgrey"]

plt.figure(figsize=(7.87,7.87))
for a,Aspect in enumerate(["mf","bp","cc"]):
    for k,Knowledge in enumerate(["NK","LK"]):
        ax=plt.subplot(3,2,2*a+k+1)

        for t,tool in enumerate(["nw","blast"]):
            fmax_list=[]
            for m,method in enumerate(method_list):
                if t==1 and method=="alnscore_1":
                    method="bitscore_1"
                infile="%s/CAFA_assessment_tool/results/%s%s_all_results.txt"%(
                    rootdir,tool,method)
                fmax=0
                if os.path.isfile(infile):
                    fp=open(infile,'r')
                    for line in fp.read().splitlines():
                        if line.startswith("%so\t%s\tfull"%(Aspect,Knowledge)):
                            fmax=float(line.split()[4])
                    fp.close()
                fmax_list.append(fmax)
                plt.text(m,fmax+0.01,("%.3f"%fmax).lstrip('0'),
                    rotation=90, fontsize=fontsize,va="bottom",
                    ha="right" if t==0 else "left")
            plt.bar(range(len(fmax_list)),fmax_list, width=0.4*(2*t-1),
                align="edge", color=color_list[t],
                label=tool)
        plt.ylabel("%s Fmax (%s targets)"%(
            Aspect.upper(),Knowledge), labelpad=0,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xticks(range(len(fmax_list)),
            [method.replace('_','') for method in method_list],
            rotation=90,fontsize=fontsize)
        plt.axis([-0.5,len(method_list)-0.5,0,0.73])
        ax.tick_params('x',length=0)
        ax.tick_params('y',length=3)
        ax.tick_params('both',pad=0)
        if a==0 and k==1:
            plt.legend(fontsize=fontsize, ncol=2, loc="upper left",
                borderpad=0.1,handlelength=1,handletextpad=0.1,
                borderaxespad=0.2,columnspacing=1)
plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
plt.savefig("Fmax_full_nw.png",dpi=300)
