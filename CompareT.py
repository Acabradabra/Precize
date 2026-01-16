#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%=================================================================================
# from unittest import case
# from IPython import get_ipython

# ip = get_ipython()
# if ip is not None:
#     ip.run_line_magic("load_ext", "autoreload")
#     ip.run_line_magic("autoreload", "2")
#%%=================================================================================
#                     Modules
#===================================================================================
import Utilities as util
(Sysa,NSysa,Arg)=util.Parseur([],0,'Arg : ')
(                             [])=Arg

from numpy import *
# import os
# import sys
# import csv
import time
import Fluent as fl

t0=time.time()
(plt,mtp)=util.Plot0()
#%%=================================================================================
#                     Parameters
#===================================================================================

Np=int(1e3)

dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/Walter-50pH2-60pJet/'
Dirs=[
dir0+'DUMP-01-Walter/' ,
dir0+'DUMP-03-A0/'     ,
dir0+'DUMP-04-Ray/'    ,
dir0+'DUMP-06-P1/' 
]

#%%=================================================================================
#                     Process
#===================================================================================
MoyT=[]
DevT=[]
for d in Dirs :
    print('=> '+d)
    dird=d+'report-meta-temperature-rfile.out'
    Dr=fl.Report_read(dird) ; Keys=list(Dr.keys()) ; print(Keys)
    MoyT.append([ mean(Dr[k][-Np:]) for k in Keys[1:] ])
    DevT.append([ std(Dr[k][-Np:]) for k in Keys[1:] ])
MoyT=array(MoyT)
DevT=array(DevT)

#%%=================================================================================
#                     Plotting
#===================================================================================

fig_V,ax_V=plt.subplots(figsize=(8,6)) #=====> Voute
Pos_V=[0.51,2.17,3.098]
labels=[ d.split('/')[-2].split('-')[-1] for d in Dirs ]
for i,l in enumerate(labels) :
    ax_V.errorbar(Pos_V,MoyT[i,:3],yerr=DevT[i,:3],marker='o',label=l)
fig_V.suptitle('Mean Voute Temperature',fontsize=20)
fig_V.legend(title='Case',fontsize=16,loc='center',bbox_to_anchor=(1.5,0.5))
ax_V.set_xlabel('X Position [m]',fontsize=16)
ax_V.set_ylabel('Temperature [K]',fontsize=16)
ax_V.set_xticks(Pos_V)
ax_V.set_xticklabels(['0.51','2.17','3.098'],fontsize=14)
# ax_V.tick_params(axis='y',labelsize=14)
# fig_V.subplots_adjust(right=2)

util.SaveFig(fig_V,'Plot/T-Voute.pdf')