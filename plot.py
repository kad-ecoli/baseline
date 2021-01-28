#!/usr/bin/env python
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import os
from os.path import dirname, basename, abspath
rootdir=dirname(abspath(__file__))

method_list=[
    "naive_1",
    "blastevalue_1",
    "blastlocalID_1",
    "blastbitscore_1",
    "blastglobalID_1",
    "blastglobalID_2",
    "blastglobalID_3",
    "blastrank_1",
    "blastfreq_1",
    "blastmetago_1",
    "blastnetgo_1",
    ]
method_list=[method for method in method_list if os.path.isfile(
    "%s/CAFA_assessment_tool/results/%s_all_results.txt"%(rootdir,method))]
fontsize=9

plt.figure(figsize=(7.87,7.87))
for a,Aspect in enumerate(["mf","bp","cc"]):
    for k,Knowledge in enumerate(["NK","LK"]):
        ax=plt.subplot(3,2,2*a+k+1)

        fmax_list=[]
        for m,method in enumerate(method_list):
            infile="%s/CAFA_assessment_tool/results/%s_all_results.txt"%(rootdir,method)
            fmax=0
            if os.path.isfile(infile):
                fp=open(infile,'r')
                for line in fp.read().splitlines():
                    if line.startswith("%so\t%s\tfull"%(Aspect,Knowledge)):
                        fmax=float(line.split()[4])
                fp.close()
            fmax_list.append(fmax)
            plt.text(m,fmax+0.01,("%.3f"%fmax).lstrip('0'),
                fontsize=fontsize,va="bottom",ha="center")
        plt.bar(range(len(fmax_list)),fmax_list,
            align="center",color="lightgrey")
        plt.ylabel("%s Fmax (%s targets)"%(
            Aspect.upper(),Knowledge), labelpad=0,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xticks(range(len(fmax_list)),
            [method.replace('_','') for method in method_list],
            rotation=90,fontsize=fontsize)
        plt.axis([-0.5,len(method_list)-0.5,0,0.7])
        ax.tick_params('x',length=0)
        ax.tick_params('y',length=3)
        ax.tick_params('both',pad=0)
plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
plt.savefig("Fmax_full.png",dpi=300)
