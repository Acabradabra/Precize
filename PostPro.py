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
(Sysa,NSysa,Arg)=util.Parseur(['Temp','Compo'],0,'Arg : ')
(                             [ TEMP , COMPO ])=Arg

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
dir0+'DUMP/'
# dir0+'DUMP-01-Walter/',
# dir0+'DUMP-03-A0/'    ,
# dir0+'DUMP-04-Ray/'   ,
# dir0+'DUMP-06-P1/'    ,
# dir0+'DUMP-07-S2S/' 
]

#===============> Param temperature
Pos_V=[0.51,2.17,3.098] # Voute thermocouples positions (m)
T_exp=array([1214,1390,1123 , 1210])+273.15 # Pilot temperatures (K) 
MoyT=[]
DevT=[]

#===============> Param composition
dil=0.5 # dilution Mf_leak/Mf_combu
Rfiles={
    'O2' : 'report-o2-rfile.out',
    'CO' : 'report-co-rfile.out',
    'CO2': 'report-co2-rfile.out',
    'H2O': 'report-h2o-rfile.out'
}
Coef= {'O2':    1e2 ,'CO':    1e6   ,'CO2': 1e2     ,'H2O':1e2      ,'N2':1e2     }
Titre={'O2':'O2 [%]','CO':'CO [ppm]','CO2':'CO2 [%]','H2O':'H2O [%]','N2':'N2 [%]'}
Keys_s=list(Rfiles.keys()) ; Nspe=len(Keys_s)
Compo_p={'O2':12.1,'CO':254,'CO2':36.25}
Compo_m={}
Compo_d={}
probe='outlet-fumes'

#%%=================================================================================
#                     Process
#===================================================================================
for d in Dirs :
    print('=> '+d)
    if TEMP : #=========================> Temperature
        f_T=d+'report-meta-temperature-rfile.out'
        Dr=fl.Report_read(f_T) ; Keys=list(Dr.keys()) ; Nl=len(Dr[Keys[0]]) ; Ns=min([Nl,Np])
        MoyT.append([ mean(Dr[k][-Ns:]) for k in Keys[1:] ])
        DevT.append([ std( Dr[k][-Ns:]) for k in Keys[1:] ])
    if COMPO : #=========================> Composition
        fig_c,ax_c=plt.subplots(figsize=(8,6),ncols=Nspe,nrows=2)
        fig_d,ax_d=plt.subplots(figsize=(8,6),ncols=Nspe)
        fig_c.suptitle('Molar fractions (%s)'%(probe),fontsize=20)
        for n,k in enumerate(Keys_s) : #=====> Wet
            f_C=d+Rfiles[k]
            Dr=fl.Report_read(f_C) ; Keys_p=list(Dr.keys()) ; Nl=len(Dr[Keys_p[0]]) ; Ns=min([Nl,Np]) #; print(Keys_p)
            key=[ kk for kk in Keys_p if probe in kk ][0]
            Compo_m[k]=mean(Dr[key][-Ns:])*Coef[k]
            Compo_d[k]=std( Dr[key][-Ns:])*Coef[k]
            ax_c[0,n].errorbar( [0] , [Compo_m[k]] , yerr=Compo_d[k], marker='o',ecolor='k',color='k')
            ax_c[0,n].set_title(Titre[k],fontsize=16)
            ax_c[0,n].set_xticks([])
        print(f"=> Wet composition : O2 {Compo_m['O2']:.2f} [%]  ,  CO {Compo_m['CO']:.0f} [ppm]  ,  CO2 {Compo_m['CO2']:.2f} [%]  ,  H2O {Compo_m['H2O']:.2f} [%]")
        C_dry=1-Compo_m['H2O']/Coef['H2O']
        X_d={ k:Compo_m[k]/C_dry for k in Keys_s[:-1] }
        for n,k in enumerate(Keys_s[:-1]) : #=====> Dry
            ax_c[1,n].errorbar( [0] , [X_d[k]] , yerr=Compo_d[k]/C_dry, marker='o',ecolor='k',color='k')
            ax_c[1,n].set_xticks([])
        print(f"=> Dry composition : O2 {X_d['O2']:.2f} [%]  ,  CO {X_d['CO']:.0f} [ppm]  ,  CO2 {X_d['CO2']:.2f} [%]")
        #====================> Dilution
        XO2_a,XN2_a=0.21,0.79
        Wa=XO2_a*32+XN2_a*28
        Wf=sum([ Compo_m[k]*fl.Mol_m[k]/Coef[k] for k in Keys_s ])
        Yf={ k:(Compo_m[k]*fl.Mol_m[k])/(Wf*Coef[k]) for k in Keys_s }    #; print('Yf :',sum([Yf[k] for k in Yf.keys()]))
        Ya={ 'O2':XO2_a*fl.Mol_m['O2']/Wa ,'N2':XN2_a*fl.Mol_m['N2']/Wa } #; print('Ya :',sum([Ya[k] for k in Ya.keys()]))
        Yd={ k:Yf[k]/(1+dil) for k in ['CO','CO2','H2O'] }
        Yd['O2']=(Yf['O2']+Ya['O2']*dil)/(1+dil)
        Yd['N2']=(         Ya['N2']*dil)/(1+dil) #; print('Yd :',sum([Yd[k] for k in Yd.keys()]))
        Wy=1/sum([ Yd[k]/fl.Mol_m[k] for k in Yd.keys() ])
        Xd_w={  k:Wy*Yd[k]/fl.Mol_m[k] for k in Yd.keys() }
        C_dry=1-Xd_w['H2O']
        Xd_d={ k:Xd_w[k]/C_dry for k in Xd_w.keys() if k!='H2O' }
        # for n,k in enumerate(Keys_s[:-1]+['N2']) : #=====> Dry + Dilution
        #     ax_d[n].plot( [0] , [Xd_d[k]*Coef[k]] , 'ok')
        #     ax_d[n].set_xticks([])
        #     # ax_c[2,n].set_title(Titre[k],fontsize=16)
        print(f"=> Dil wet composition : O2 {Xd_w['O2']*Coef['O2']:.2f} [%]  ,  CO {Xd_w['CO']*Coef['CO']:.0f} [ppm]  ,  CO2 {Xd_w['CO2']*Coef['CO2']:.2f} [%]  ,  N2 {Xd_w['N2']*Coef['N2']:.2f} [%]  ,  H2O {Xd_w['H2O']*Coef['H2O']:.2f} [%]")
        print(f"=> Dil dry composition : O2 {Xd_d['O2']*Coef['O2']:.2f} [%]  ,  CO {Xd_d['CO']*Coef['CO']:.0f} [ppm]  ,  CO2 {Xd_d['CO2']*Coef['CO2']:.2f} [%]  ,  N2 {Xd_d['N2']*Coef['N2']:.2f} [%]")
        ax_c[1,-1].axis('off')
        ax_c[0,0].set_ylabel('Wet composition',fontsize=16)
        ax_c[1,0].set_ylabel('Dry composition',fontsize=16)
        util.SaveFig( fig_c,d+'Plot/Compo.pdf' )

MoyT=array(MoyT)
DevT=array(DevT)

#%%=================================================================================
#                     Plotting Temperature
#===================================================================================
if TEMP :
    fig_V,ax_V=plt.subplots(figsize=(8,6)) #=====> Voute
    fig_C,ax_C=plt.subplots(figsize=(8,6)) #=====> Chimney
    labels=[ d.split('/')[-2].split('-')[-1] for d in Dirs ] ; Nl=len(labels)
    for i,l in enumerate(labels) :
        ax_C.plot([i,i], [T_exp[ -1],MoyT[i,-1]],'k:',alpha=0.3)
        Er=100*abs(MoyT[i,-1]-T_exp[ -1])/T_exp[ -1]
        ax_C.text(i+0.1,T_exp[ -1]+50,f'{Er:.1f} %',fontsize=12,color='k')
        ax_V.errorbar(Pos_V,MoyT[i,:3],yerr=DevT[i,:3],marker='o',label=l)
        ax_C.errorbar(i    ,MoyT[i,-1],yerr=DevT[i,-1],marker='o',label=l)
    ax_V.plot(Pos_V ,   T_exp[:-1] ,'k-o',label='Pilote')
    ax_C.plot(   Nl ,   T_exp[ -1] ,'k-o',label='Pilote')
    ax_C.plot([0,Nl],2*[T_exp[ -1]],'k--')

    fig_V.suptitle('Mean Voute Temperature',fontsize=20)
    fig_C.suptitle('Mean Chimney Temperature',fontsize=20)
    fig_V.legend(title='Case',fontsize=16,loc='lower left',bbox_to_anchor=(0.35,0.12))
    ax_V.set_xlabel('X Position [m]'  ,fontsize=16)
    ax_V.set_ylabel('Temperature [K]' ,fontsize=16)
    ax_C.set_xlabel('Radiation models',fontsize=16)
    ax_C.set_ylabel('Temperature [K]' ,fontsize=16)
    ax_V.set_xticks(Pos_V)
    ax_C.set_xticks(range(Nl+1))
    ax_V.set_xticklabels(['0.51','2.17','3.098'],fontsize=14)
    ax_C.set_xticklabels(labels+['Pilote'],fontsize=14)
    fig_V.subplots_adjust(right=0.5)

    util.SaveFig(fig_V,'Plot/T-Voute.pdf')
    util.SaveFig(fig_C,'Plot/T-Chemine.pdf')