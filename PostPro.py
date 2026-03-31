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
(Sysa,NSysa,Arg)=util.Parseur(['Temp','Compo','Report','Visu','Prof','Voute'],0,'Arg : ')
(                             [ TEMP , COMPO , REPORT , VISU , PROF , VOUTE ])=Arg

from numpy import *
import os
# import sys
# import csv
import time
import Fluent as fl

import ParamPilotV2 as pp

t0=time.time()
(plt,mtp)=util.Plot0()
#%%=================================================================================
#                     Parameters
#===================================================================================

Np=int(1e3)

# dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/Walter-50pH2-60pJet/'
# dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M00/'
dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M01/'
# dir0='/scratch/ZEUS/FLUENT/PRECIZE/RUN-M01/'
Dirs=[
# dir0+'DUMP/'
dir0+'DUMP-16-OD2-Walls/'
]

#===============> Param Visualisation
# Vars=['Vel','k','T','o2','h2','ch4','co2','h2o','co']
cmesh=1e3 # coef coordinate (0 : no ticks)
Keys=['Tf']
Surfs  ={'Tf':'Data-F-x0.dat'  }
Planes ={'Tf':'yz'             }
Vars   ={'Tf':'temperature'    }
Titres ={'Tf':'Temperature [K]'}
RX     ={'Tf':[]               }
RY     ={'Tf':[-0.2,2]         }
Ticks  ={'Tf':[]               }
BD_Vars={'Tf':[]               }
Fsizes ={'Tf':(15,7)           }
Cmaps  ={'Tf':'inferno'        }
Names  ={'Tf':'Temperature.png'}
Params ={'Tf':[]               }
#===============> Param profile
fprof='Data-F-chem.dat'
fvout='Data-TOP.dat'

#===============> Param temperature
Pos_V=[0.51,2.17,3.098] # Voute thermocouples positions (m)
# T_exp=array([1214,1390,1123 , 1210])+273.15 # Pilot temperatures (K) 
# T_exp=array([1159.1,1351,1110.8 , 1170.2])+273.15 # Pilot temperatures (K) 
T_exp=array([1150,1400,1150 , 1350])+273.15 # Pilot temperatures (K) 
MoyT=[]
DevT=[]

#===============> Param composition
# dil=0.45 # dilution Mf_leak/Mf_tot
# dil=0.20 # dilution Mf_leak/Mf_tot
dil=0 # dilution Mf_leak/Mf_tot
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
Compo_p['N2']=100-(Compo_p['O2']+Compo_p['CO']*1e-4+Compo_p['CO2'])
Compo_m={}
Compo_d={}
# probe='outlet-fumes'
probe='pd'
# probe='sample_d'
# probe='sample_c'

#===============> Inlet properties
Tu=300 # Ambient temperature (K)

#%%=================================================================================
#                     Process
#===================================================================================
if PROF :
    figP,axP=plt.subplots(figsize=(10,6))
T_dil=[]
for d in Dirs :
    print('=> '+d) ; dirp=d+'PLOT/'
    if TEMP : #===========================================================================> Temperature
        f0=d+'report'
        Dr=fl.Report_read(f0+'-temperature-rfile.out') ; Keys=list(Dr.keys()) ; Nl=len(Dr[Keys[0]]) ; Ns=min([Nl,Np]) ; KeysT=Keys
        Id_v=[ KeysT.index(k) for k in ["p3","p2","p1"] ]
        I_out=KeysT.index('zc')
        MoyT.append([ mean(Dr[k][-Ns:]) for k in KeysT ]) ; Tb=float(MoyT[-1][I_out]) ; print(f'=> Mean burned temperature : {Tb:.0f} [K]')
        DevT.append([ std( Dr[k][-Ns:]) for k in KeysT ])
        Dq=fl.Report_read(d+'report-heat-release-rfile.out') ; Keys=list(Dq.keys())
        Hr=mean(Dq[Keys[1]][-Ns:]) ; print(f'=> Mean heat release : {Hr*1e-3:.0f} [kW]')
        Dm=fl.Report_read(d+'report-mass-balance-rfile.out') ; Keys=list(Dm.keys())
        (Mf_f,Mf_o,Mf_b,Mf_l,Mf_s,Mb)=fl.Mf_sep2(Dm,Keys) ; Ns=min([len(Mb),Np]) ; Md=mean(Mf_b[-Ns:])
        Cp=-Hr/(Md*(Tb-Tu)) ; print(f'=> Estimated Cp : {Cp:.2f} J/Kg/K')
        Tb2=Tu-Hr/((1+dil)*Md*Cp) ; print(f'=> Estimated Tb (with dilution {dil*1e2:.0f} %) : {Tb2:.2f} [K]')
        T_dil.append(Tb2)
    if COMPO : #===========================================================================> Composition
        fig_c,ax_c=plt.subplots(figsize=(8,6),ncols=Nspe,nrows=2)
        fig_d,ax_d=plt.subplots(figsize=(8,6),ncols=Nspe)
        fig_c.suptitle('Molar fractions (%s)'%(probe),fontsize=20)
        fig_d.suptitle(f'Dry molar fraction : {dil*1e2:.0f} % diluted',fontsize=20)
        for n,k in enumerate(Keys_s) : #=====> Wet
            f_C=d+Rfiles[k]
            if os.path.exists(f_C) :
                Dr=fl.Report_read(f_C) ; Keys_p=list(Dr.keys()) ; Nl=len(Dr[Keys_p[0]]) ; Ns=min([Nl,Np]) #; print(Keys_p)
                key=[ kk for kk in Keys_p if probe in kk ][0]
                Compo_m[k]=mean(Dr[key][-Ns:])*Coef[k]
                Compo_d[k]=std( Dr[key][-Ns:])*Coef[k]
            else :
                print(util.Col('r',f'=> No report file found for {k} !')) ; Compo_m[k]=0 ; Compo_d[k]=0
            ax_c[0,n].errorbar( [0] , [Compo_m[k]] , yerr=Compo_d[k], marker='o',ecolor='k',color='k')
            ax_c[0,n].set_title(Titre[k],fontsize=16)
            ax_c[0,n].set_xticks([])
        print(f"=> Wet composition : O2 {Compo_m['O2']:.2f} [%]  ,  CO {Compo_m['CO']:.0f} [ppm]  ,  CO2 {Compo_m['CO2']:.2f} [%]  ,  H2O {Compo_m['H2O']:.2f} [%]")
        C_dry=1-Compo_m['H2O']/Coef['H2O']
        X_d={ k:Compo_m[k]/C_dry for k in Keys_s[:-1] }
        for n,k in enumerate(Keys_s[:-1]) : #=====> Dry
            ax_c[1,n].errorbar( [0] , [X_d[k]] , yerr=Compo_d[k]/C_dry, marker='o',ecolor='k',color='k')
            ax_c[1,n].set_xticks([])
        print(util.Col('b',f"=> Dry composition : O2 {X_d['O2']:.2f} [%]  ,  CO {X_d['CO']:.0f} [ppm]  ,  CO2 {X_d['CO2']:.2f} [%]"))
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
        for n,k in enumerate(Keys_s[:-1]+['N2']) : #=====> Dry + Dilution
            E=100*abs(Compo_p[k]-Xd_d[k]*Coef[k])/Compo_p[k]
            ax_d[n].plot(2*[0] , [Compo_p[k],Xd_d[k]*Coef[k]] , ':k')
            ax_d[n].plot(  [0] , [   Xd_d[k]*Coef[k]] , 'or')
            ax_d[n].plot(  [0] , [Compo_p[k]        ] , 'ok')
            ax_d[n].set_xticks([])
            ax_d[n].set_title(Titre[k],fontsize=16)
            ax_d[n].set_xlabel(f'{E:.1f} % error',fontsize=12)
        print(f"=> Dil wet composition : O2 {Xd_w['O2']*Coef['O2']:.2f} [%]  ,  CO {Xd_w['CO']*Coef['CO']:.0f} [ppm]  ,  CO2 {Xd_w['CO2']*Coef['CO2']:.2f} [%]  ,  N2 {Xd_w['N2']*Coef['N2']:.2f} [%]  ,  H2O {Xd_w['H2O']*Coef['H2O']:.2f} [%]")
        print(f"=> Dil dry composition : O2 {Xd_d['O2']*Coef['O2']:.2f} [%]  ,  CO {Xd_d['CO']*Coef['CO']:.0f} [ppm]  ,  CO2 {Xd_d['CO2']*Coef['CO2']:.2f} [%]  ,  N2 {Xd_d['N2']*Coef['N2']:.2f} [%]")
        ax_c[1,-1].axis('off')
        ax_c[0,0].set_ylabel('Wet composition',fontsize=16)
        ax_c[1,0].set_ylabel('Dry composition',fontsize=16)
        util.SaveFig( fig_c,d+'Plot/Compo-%s.pdf'%(probe) )
        util.SaveFig( fig_d,d+'Plot/Compo-%s_diluted.pdf'%(probe) )
    if REPORT : #===========================================================================> Report
        for rf in os.popen('ls %s/report-*-rfile.out'%(d)).read().splitlines() :
            r_name=rf.split('/')[-1][7:-10]
            Dr=fl.Report_read(rf) ; Keys=list(Dr.keys()) #; print(Keys)
            if '(' in Keys[1] : Labels=['Iteration']+[ k.split('(')[1][:-1] for k in Keys[1:] if '(' in k ]
            else              : Labels=Keys
            print('=> \033[31m%s\033[0m : '%(r_name) , Labels )
            figr,axr=plt.subplots(figsize=(10,7)) #; bxr=axr.twinx()
            if r_name == 'mass-balance' :
                (Mf_f,Mf_o,Mf_b,Mf_l,Mf_s,Mb)=fl.Mf_sep2(Dr,Keys) ; Ns=min([len(Mb),Np])
                Ma_f=mean(Mf_f[-Ns:])
                Ma_o=mean(Mf_o[-Ns:])
                Ma_b=mean(Mf_b[-Ns:])
                Ma_l=mean(Mf_l[-Ns:])
                Ma_s=mean(Mf_s[-Ns:])
                Mbal=mean(  Mb[-Ns:])
                Ma_t=Ma_f+Ma_o+Ma_s
                Q_l=3600*Ma_l/fl.Rho0_air ; R_l=100*Ma_l/Ma_t
                print(f'=> fuel : {Ma_f*1e3:.3f} g/s  ,  oxid : {Ma_o*1e3:.3f} g/s  ,  out : {Ma_b*1e3:.3f} g/s  ,   leak : {Ma_l*1e3:.3f} g/s  ,  slope : {Ma_s*1e3:.3f} g/s  ,  balance : {Mbal*1e3:.3f} g/s')
                print(f'=> Leak : {Q_l:.3f} [Nm3/h]  ,  {R_l:.2f} % M tot')
                print(f'=> Balance : {-100*Mbal/Ma_t:.2f} [%] O2 + fuel + slope')
                axr.plot( Dr['Iteration'],1e3*Mf_f  ,label='inlet-fuel'  )
                axr.plot( Dr['Iteration'],1e3*Mf_o  ,label='inlet-oxid'  )
                axr.plot( Dr['Iteration'],1e3*Mf_b  ,label='outlet'      )
                axr.plot( Dr['Iteration'],1e3*Mf_l  ,label='leak'        )
                axr.plot( Dr['Iteration'],1e3*Mf_s  ,label='slope-zone'  )
                # bxr.plot( Dr['Iteration'],1e3*Mb,'k',label='mass-balance')
                axr.plot( Dr['Iteration'],1e3*Mb,'k',label='mass-balance')
                axr.set_ylabel('Mass flow rate [g/s]',fontsize=25)
                # figr.legend(fontsize=15,loc='center',bbox_to_anchor=(0.7,0.3))
                # axr.legend(fontsize=15,loc='best')
                axr.legend(fontsize=15,loc='lower right')
                #======> Detail
                Ndet=250
                figd,axd=plt.subplots(ncols=2,figsize=(12,6))
                axd[0].set_title('Leaks',fontsize=25)
                axd[0].plot( Dr['Iteration'][-Ndet:],1e3*Dr['f-bec-in'][-Ndet:],label='bec'    )
                axd[0].plot( Dr['Iteration'][-Ndet:],1e3*Dr['pa-in'   ][-Ndet:],label='pyro-a' )
                axd[0].plot( Dr['Iteration'][-Ndet:],1e3*Dr['pc-in'   ][-Ndet:],label='pyro-c' )
                axd[0].plot( Dr['Iteration'][-Ndet:],1e3*Dr['pd-in'   ][-Ndet:],label='pyro-d' )
                axd[0].plot( Dr['Iteration'][-Ndet:],1e3*Dr['v-in'    ][-Ndet:],label='visse'  )
                axd[1].set_title('Outlet',fontsize=25)
                axd[1].plot( Dr['Iteration'][-Ndet:],1e3*Dr['f-out'   ][-Ndet:],label='outlet')
                axd[1].plot( Dr['Iteration'][-Ndet:],1e3*Dr['f-jeu'   ][-Ndet:],label='jeu'   )
                axd[1].plot( Dr['Iteration'][-Ndet:],1e3*Dr['f-jupe'  ][-Ndet:],label='jupe'  )
                axd[1].plot( Dr['Iteration'][-Ndet:],1e3*Dr['f-dilute'][-Ndet:],label='dilute')
                axd[0].set_ylabel('Mass flow rate [g/s]',fontsize=25)
                axd[0].legend(fontsize=15)
                axd[1].legend(fontsize=15)
                util.SaveFig(figd,dirp+'MassFlow-Details.pdf')
            elif r_name == 'heat-release' :
                for n,k in enumerate(Keys[1:]) : 
                    if 'fluides-fluide' in k :
                        axr.plot( Dr['Iteration'],Dr[k],label=Labels[n+1] )
            elif r_name == 'heat-loss' :
                Nskip=20
                Ns=min([len(Dr['Iteration']),Np])
                figW,axW=plt.subplots(figsize=(10,7))
                figW.suptitle('Wall heat losses',fontsize=20)
                axr.set_ylabel('Heat loss [kW]',fontsize=25)
                axW.set_ylabel('Heat loss [kW]',fontsize=25)
                for n,k in enumerate(Keys[1:]) : 
                        if k not in ['f-ts','f-tt','f-vs'] :
                            axr.plot( Dr['Iteration'][Nskip:],Dr[k][Nskip:]*1e-3,label=Labels[n+1] )
                        if k in ['f-tc'] : 
                            hl_m=mean(Dr[k][-Ns:]) # 
                            print(f'=> {k} : {hl_m*1e-3:.0f} [kW]  ,  {100*hl_m/(pp.Pow*1e3):.2f} % Pow')
                        if k in ['f-top','f-side','f-front','f-back'] :
                            axW.plot( Dr['Iteration'][Nskip:],Dr[k][Nskip:]*1e-3,label=Labels[n+1] )
                if len(Keys)>2 : 
                    axr.legend(fontsize=15)
                    axW.legend(fontsize=15)
                (hl_walls,hl_talus,hl_wb)=fl.Hl_sep(Dr,Ns)
                print(f'=> Total wall heat loss : {hl_walls*1e-3:.0f} [kW]  ,  {100*hl_walls/(pp.Pow*1e3):.2f} % Pow')
                print(f'=> Total talu heat loss : {hl_talus*1e-3:.0f} [kW]  ,  {100*hl_talus/(pp.Pow*1e3):.2f} % Pow')
                print(f'=> Water boxe heat loss : {hl_wb   *1e-3:.0f} [kW]  ,  {100*hl_wb   /(pp.Pow*1e3):.2f} % Pow')
                util.SaveFig(figW,dirp+'HeatLoss-Details.pdf')
            elif r_name == 'temperature' :
                figE,axE=plt.subplots(figsize=(10,7))
                figE.suptitle('External wall temperature',fontsize=20)
                figr.suptitle('Temperature [K]',fontsize=20)
                for n,k in enumerate(Keys[1:]) : 
                    if '-2:external' not in k : axr.plot( Dr['Iteration'],Dr[k],label=Labels[n+1] )
                    else                      : axE.plot( Dr['Iteration'],Dr[k],label=Labels[n+1] )
                if len(Keys)>2 : 
                    axr.legend(fontsize=15)
                    axE.legend(fontsize=15)
                util.SaveFig(figE,dirp+'T-External.pdf')
            else :
                for n,k in enumerate(Keys[1:]) : axr.plot( Dr['Iteration'],Dr[k],label=Labels[n+1] )
                if len(Keys)>2 : axr.legend(fontsize=15) #,loc='center',bbox_to_anchor=(0.7,0.3))
            # if len(Keys)>2 : figr.legend(fontsize=15,loc='center',bbox_to_anchor=(0.7,0.3))
            if r_name in ['heat-release','mass-balance','shell-loss','co2'] :
                axr.ticklabel_format( axis='y' , scilimits=(-3,3) )
            util.SaveFig(figr,dirp+'Report-%s.pdf'%(r_name))
    if VOUTE : #===========================================================================> Voute
        Data=fl.ReadSurfD(d+'DATA/'+fvout)
    if PROF : #===========================================================================> profile
        Data=fl.ReadSurfD(d+'DATA/'+fprof)
    if VISU : #===========================================================================> Visualisation
        F_int={}
        for k in Vars :
            F_int[k]=fl.Visu(d+'DATA/'+Surfs[k],Planes[k],Vars[k],Titres[k],RX[k],RY[k],Ticks[k],cmesh,BD_Vars[k],Fsizes[k],Cmaps[k],dirp+Names[k],Params[k])

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
        Tchem=MoyT[i,I_out]
        ax_C.plot([i,i],[T_exp[-1],Tchem],'k:',alpha=0.3)
        ax_C.plot([i  ], T_dil[i],'ko')
        Er=100*abs(Tchem-T_exp[ -1])/T_exp[ -1]
        ax_C.text(i+0.1,T_exp[ -1]+50,f'{Er:.1f} %',fontsize=12,color='k')
        ax_V.errorbar(Pos_V,MoyT[i,Id_v],yerr=DevT[i,Id_v ],marker='o',label=l)
        ax_C.errorbar(i    ,Tchem       ,yerr=DevT[i,I_out],marker='o',label=l)
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

    if len(Dirs)>1 : d_compa=dir0
    else           : d_compa=Dirs[0]
    util.SaveFig(fig_V,d_compa+'PLOT/T-Voute.pdf')
    util.SaveFig(fig_C,d_compa+'PLOT/T-Chemine.pdf')