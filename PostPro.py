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
(Sysa,NSysa,Arg)=util.Parseur(['Temp','Compo','Report','Visu','Prof','Voute','Dump','All'],0,'Arg : ')
(                             [ TEMP , COMPO , REPORT , VISU , PROF , VOUTE , DUMP , ALL ])=Arg
if ALL : TEMP,COMPO,REPORT,VISU,PROF,VOUTE=True,True,True,True,True,True

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

#===============> Directories
# dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/Walter-50pH2-60pJet/'
# dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M00/'
# dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M01/'
# dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M02/'
dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M03/'
# dir0='/scratch/ZEUS/FLUENT/PRECIZE/RUN-M01/'
Dirs=[
# 'DUMP-00-Init/'
# 'DUMP-01-Jeu/'
'DUMP-01-Pstat/'
# 'DUMP-18-OD2-Bath/'
]
if DUMP : Dirs=[ 'DUMP/' ]
if len(Dirs)>1 : d_compa=dir0
else           : d_compa=dir0+Dirs[0]

#===============> Param Visualisation
cmesh=0 # coef coordinate (0 : no ticks)
# cmesh=1e2 # coef coordinate (0 : no ticks)
# Vars=['Tf','YCO2f','YH2Of','Ptf','Psf','Velf']
Vars=['Vcol','Tcol']
Param_visu={
'Tcol' :['Data-F-yc.dat','xz','temperature'       ,'Temperature [K]'     ,[-0.3,0.3],[2,4]   ,[],cmesh,[]         ,(5,9) ,'inferno','Col-Temperature.png' ,['ISO',[1800]]],
'Vcol' :['Data-F-yc.dat','xz','velocity-magnitude','Velocity [m/s]'      ,[-0.3,0.3],[2,3]   ,[],cmesh,[]         ,(7,9) ,'cividis','Col-Velocity.png'    ,['QUIV',[2e-2,2e-2,300]]],
'YCO2f':['Data-F-x0.dat','yz','co2'               ,'$Y_{CO_2}$ [-]'      ,[]        ,[-0.2,2],[],cmesh,[]         ,(15,7),'viridis','Four-YCO2.png'       ,[]            ],
'YH2Of':['Data-F-x0.dat','yz','h2o'               ,'$Y_{H_2O}$ [-]'      ,[]        ,[-0.2,2],[],cmesh,[]         ,(15,7),'viridis','Four-YH2O.png'       ,[]            ],
'Tf'   :['Data-F-x0.dat','yz','temperature'       ,'Temperature [K]'     ,[]        ,[-0.2,2],[],cmesh,[ 290,2500],(15,7),'inferno','Four-Temperature.png',['ISO',[2000,2500]]],
'Ptf'  :['Data-F-x0.dat','yz','total-pressure'    ,'Total pressure [Pa]' ,[]        ,[-0.2,2],[],cmesh,[]         ,(15,7),'cividis','Four-Pt.png'         ,['ISO',[0]]        ],
'Psf'  :['Data-F-x0.dat','yz','pressure'          ,'Static pressure [Pa]',[]        ,[-0.2,2],[],cmesh,[]         ,(15,7),'cividis','Four-Ps.png'         ,['ISO',[0]]        ],
'Velf' :['Data-F-x0.dat','yz','velocity-magnitude','Velocity [m/s]'      ,[]        ,[-0.2,2],[],cmesh,[]         ,(15,7),'cividis','Four-Velocity.png'   ,['QUIV',[1e-1,1e-1,300]]],
'Ttop' :['Data-TOP.dat' ,'yx','temperature'       ,'Temperature [K]'     ,[]        ,[]      ,[],cmesh,[]         ,(15,7),'inferno','Top-Temperature.png' ,['INTERP']    ],
'Htop' :['Data-TOP.dat' ,'yx','heat-flux'         ,'Heat flux [W/m2]'    ,[]        ,[]      ,[],cmesh,[]         ,(15,7),'inferno','Top-HeatFlux.png'    ,['INTERP']    ],
}
if ALL : Vars=list(Param_visu.keys())
if dir0=='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M01/' :
    Param_visu['Ptf' ][8]=[-30,20]
    Param_visu['Psf' ][8]=[-40,1e-9]
    Param_visu['Velf'][8]=[  0,25]
    Param_visu['Ttop'][8]=[1740,1790]
    Param_visu['Htop'][8]=[-8e3,-2e3]
    Param_visu['Tcol'][8]=[299,1800]
    Param_visu['Tcol'][-1][-1]=[1000]
elif dir0=='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M02/' :
    Param_visu['Velf'][8]=[   0, 25]
    Param_visu['Ttop'][8]=[]
    Param_visu['Htop'][8]=[]
    if Dirs[0][5:7]=='00' :
        Param_visu['Psf' ][8]=[  70, 95]
        Param_visu['Ptf' ][8]=[  70,110]
    if Dirs[0][5:7]=='01' :
        Param_visu['Psf' ][8]=[-5,25]
        Param_visu['Ptf' ][8]=[-5,50]

#===============> Param Voute
dh=25e-3 # [m] Height of thermocouples
lkk=2.8  # [W/m K] Radex conductivity

#===============> Param profile
fprof='Data-F-chem.dat'
Z_four=[1.344,2.015,2.117,2.572,2.768 , 12.5]
Ps_top=-pp.Rho_air*9.81*Z_four[-1]
Param_prof={
    'T'   :['temperature'   ,'Temperature [K]'     ,[ 600,2000],'Profile-T.pdf' ],
    'Pt'  :['total-pressure','Total pressure [Pa]' ,[-150, 100],'Profile-Pt.pdf'],
    'Ps'  :['pressure'      ,'Static pressure [Pa]',[-200, 100],'Profile-Ps.pdf'],
    'YCO2':['co2'           ,'$Y_{CO_2}$ [-]'      ,[   0, 0.5],'Profile-YCO2.pdf'],
}

#===============> Param temperature
Pos_V=[0.51,2.17,3.098] # Voute thermocouples positions (m)
# Pos_V=[0.51,1.377,3.098] # Voute thermocouples positions (m)
# T_exp=array([1214,1390,1123 , 1210])+273.15 # Pilot temperatures (K) 
# T_exp=array([1159.1,1351,1110.8 , 1170.2])+273.15 # Pilot temperatures (K) 
T_exp=array([ 1150,1400,1150 , 1350 ])+273.15 # Pilot temperatures (K) 
MoyT=[]
DevT=[]

#===============> Param Heat losses
Nskip=100

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
Keys_s=list(Rfiles.keys())+['N2'] ; Nspe=len(Keys_s)+1
Keys_w=[ 'O2' , 'CO' , 'CO2' , 'N2' , 'H2O' ]
Keys_d=[ k for k in Keys_s if k!='H2O' ]
# Compo_p={'O2':12.1,'CO':254,'CO2':36.25}
# Compo_p={'O2':12,'CO':774,'CO2':38} # 100% H2 L255
Compo_p={'O2':6,'CO':1894,'CO2':41} # 100% H2 L254
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
if VOUTE :
    if VISU==False : VISU=True
    if 'Ttop' not in Vars : Vars.append('Ttop')
    if 'Htop' not in Vars : Vars.append('Htop')
    figV,axV=plt.subplots(figsize=(10,6)) #; bxV=axV.twinx()
    figV.suptitle('Voute Temperature',fontsize=20)
    axV.set_xlabel('Y position [m]',fontsize=16)
    axV.set_ylabel('Temperature [K]',fontsize=16)
    axV.plot(Pos_V,T_exp[:3],'ok',label='Pilot')
#===================================================================================
if PROF :
    FIG_P,AX_P={},{}
    for k in Param_prof.keys() :
        figP,axP=plt.subplots(figsize=(10,6))
        figP.suptitle('Vertical profile through chimney',fontsize=20)
        axP.set_xlabel('Z position [m]',fontsize=16)
        axP.set_ylabel(Param_prof[k][1],fontsize=16)
        axP.set_ylim(Param_prof[k][2])
        for z in Z_four[:-1] : axP.plot( 2*[z],Param_prof[k][2],':k')
        FIG_P[k]=figP ; AX_P[k]=axP
    AX_P['Pt'].plot([0,Z_four[-1]],[0,Ps_top],'k' )
    AX_P['Ps'].plot([0,Z_four[-1]],[0,Ps_top],'k' )
#===================================================================================
T_dil=[]
for d0 in Dirs :
    d=dir0+d0
    print('=> '+d) ; dirp=d+'PLOT/'
    if TEMP   : #===========================================================================> Temperature
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
    if COMPO  : #===========================================================================> Composition
        fig_c,ax_c=plt.subplots(figsize=(10,6),ncols=Nspe-1,nrows=2)
        fig_d,ax_d=plt.subplots(figsize=( 8,6),ncols=Nspe-2)
        fig_c.suptitle('Molar fractions (%s)'%(probe),fontsize=20)
        # fig_d.suptitle(f'Dry molar fraction ({'\033[31m'}Simu{'\033[0m'}Pilot) : {dil*1e2:.0f} % diluted',fontsize=20)
        fig_d.suptitle(f'Dry molar fraction : {dil*1e2:.0f} % diluted',fontsize=20)
        for n,k in enumerate(Rfiles) : #=====> Wet
            f_C=d+Rfiles[k]
            if os.path.exists(f_C) :
                Dr=fl.Report_read(f_C) ; Keys_p=list(Dr.keys()) ; Nl=len(Dr[Keys_p[0]]) ; Ns=min([Nl,Np]) #; print(Keys_p)
                key=[ kk for kk in Keys_p if probe in kk ][0]
                Compo_m[k]=mean(Dr[key][-Ns:])*Coef[k]
                Compo_d[k]=std( Dr[key][-Ns:])*Coef[k]
            else :
                print(util.Col('r',f'=> No report file found for {k} !')) ; Compo_m[k]=0 ; Compo_d[k]=0
        Compo_m['N2']=(1-(Compo_m['O2']/Coef['O2']+Compo_m['CO']/Coef['CO']+Compo_m['CO2']/Coef['CO2']+Compo_m['H2O']/Coef['H2O']))*Coef['N2']
        Compo_d['N2']=mean([ Compo_d[k]/Coef[k] for k in Keys_s[:-1] ])*Coef['N2']
        for n,k in enumerate(Keys_w) :
            ax_c[0,n].errorbar( [0] , [Compo_m[k]] , yerr=Compo_d[k], marker='o',ecolor='k',color='k')
            ax_c[0,n].set_title(Titre[k],fontsize=16)
            ax_c[0,n].set_xticks([])
        print(f"=> Wet composition : O2 {Compo_m['O2']:.2f} [%]  ,  CO {Compo_m['CO']:.0f} [ppm]  ,  CO2 {Compo_m['CO2']:.2f} [%]  ,  H2O {Compo_m['H2O']:.2f} [%]  ,  N2 {Compo_m['N2']:.2f} [%]")
        C_dry=1-Compo_m['H2O']/Coef['H2O']
        X_d={ k:Compo_m[k]/C_dry for k in Keys_s }
        for n,k in enumerate(Keys_d) : #=====> Dry
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
        Yd['N2']=(Yf['N2']+Ya['N2']*dil)/(1+dil) #; print('Yd :',sum([Yd[k] for k in Yd.keys()]))
        Wy=1/sum([ Yd[k]/fl.Mol_m[k] for k in Yd.keys() ])
        Xd_w={  k:Wy*Yd[k]/fl.Mol_m[k] for k in Yd.keys() }
        C_dry=1-Xd_w['H2O']
        Xd_d={ k:Xd_w[k]/C_dry for k in Xd_w.keys() if k!='H2O' }
        for n,k in enumerate(Keys_d) : #=====> Dry + Dilution
            E=100*abs(Compo_p[k]-Xd_d[k]*Coef[k])/Compo_p[k]
            ax_d[n].plot(2*[0] , [Compo_p[k],Xd_d[k]*Coef[k]] , ':k')
            ax_d[n].plot(  [0] , [   Xd_d[k]*Coef[k]] , 'or')
            ax_d[n].plot(  [0] , [Compo_p[k]        ] , 'ok')
            ax_d[n].set_xticks([])
            ax_d[n].set_yticks([ round(x, 1) for x in sorted([Xd_d[k]*Coef[k],Compo_p[k]]) ])
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
                if 'f-bec-in' in Dr : axd[0].plot( Dr['Iteration'][-Ndet:],1e3*Dr['f-bec-in'][-Ndet:],label='bec'    )
                if 'pa-in'    in Dr : axd[0].plot( Dr['Iteration'][-Ndet:],1e3*Dr['pa-in'   ][-Ndet:],label='pyro-a' )
                if 'pc-in'    in Dr : axd[0].plot( Dr['Iteration'][-Ndet:],1e3*Dr['pc-in'   ][-Ndet:],label='pyro-c' )
                if 'pd-in'    in Dr : axd[0].plot( Dr['Iteration'][-Ndet:],1e3*Dr['pd-in'   ][-Ndet:],label='pyro-d' )
                if 'v-in'     in Dr : axd[0].plot( Dr['Iteration'][-Ndet:],1e3*Dr['v-in'    ][-Ndet:],label='visse'  )
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
                            # print(f'=> {k} : {hl_m*1e-3:.0f} [kW]  ,  {100*hl_m/(pp.Pow*1e3):.2f} % Pow')
                        if k in ['f-top','f-side','f-front','f-back'] :
                            axW.plot( Dr['Iteration'][Nskip:],Dr[k][Nskip:]*1e-3,label=Labels[n+1] )
                if len(Keys)>2 : 
                    axr.legend(fontsize=15)
                    axW.legend(fontsize=15)
                (hl_walls,hl_talus,hl_wb,hl_bath)=fl.Hl_sep(Dr,Ns,pp.Pow,Verbose=2)
                util.SaveFig(figW,dirp+'HeatLoss-Details.pdf')
            elif r_name == 'temperature' :
                figE,axE=plt.subplots(figsize=(10,7))
                figE.suptitle('External wall temperature [K]',fontsize=20)
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
    if PROF   : #===========================================================================> profile
        Data=fl.ReadSurfD(d+'DATA/'+fprof)
        for k in Param_prof : AX_P[k].plot(Data['z-coordinate'],Data[Param_prof[k][0]])
    if VISU   : #===========================================================================> Visualisation
        F_int={}
        for k in Vars : 
            Param_visu[k][0 ]=d+'DATA/'+Param_visu[k][0 ]
            Param_visu[k][11]=dirp     +Param_visu[k][11]
            F_int[k]=fl.Visu(*Param_visu[k])
    if VOUTE  : #===========================================================================> Voute
        Vy=linspace(0.4,3.7,int(1e3))
        Vx=0*Vy
        T_int=F_int['Ttop'](Vy,Vx)
        H_int=F_int['Htop'](Vy,Vx)
        Vy2=Vy-Vy[0]
        T_w1  =T_int+(dh+1e-3)*H_int/lkk
        T_wall=T_int+ dh      *H_int/lkk
        T_w2  =T_int+(dh-1e-3)*H_int/lkk
        # axV.plot(Vy,T_int ,'k')
        axV.plot(Vy,T_wall,'b')
        axV.plot(Vy,T_w1  ,':b')
        axV.plot(Vy,T_w2  ,':b')
        # bxV.plot(Vy,H_int,'r')
        util.SaveFig(figV,dirp+'T-Voute-Int.pdf')

if PROF : 
    for k in Param_prof.keys() : util.SaveFig(FIG_P[k],d_compa+Param_prof[k][-1])

MoyT=array(MoyT)
DevT=array(DevT)

#%%=================================================================================
#                     Plotting Temperature
#===================================================================================
if TEMP :
    fig_V,ax_V=plt.subplots(figsize=(8,6)) #=====> Voute
    fig_C,ax_C=plt.subplots(figsize=(8,6)) #=====> Chimney
    labels=[ d.split('/')[-2].split('-')[-1] for d in Dirs ] ; Nl=len(labels)
    if Nl==1 : labels=['Simu']
    for i,l in enumerate(labels) :
        Tchem=MoyT[i,I_out]
        ax_C.plot([i,i],[T_exp[-1],Tchem],'k:',alpha=0.3)
        ax_C.plot([i  ], T_dil[i],'ko')
        Er=100*abs(Tchem-T_exp[ -1])/T_exp[ -1]
        ax_C.text(i+0.1,0.5*(T_exp[ -1]+Tchem),f'{Er:.1f} %',fontsize=12,color='k')
        ax_V.errorbar(Pos_V,MoyT[i,Id_v],yerr=DevT[i,Id_v ],marker='o',label=l)
        ax_C.errorbar(i    ,Tchem       ,yerr=DevT[i,I_out],marker='o',label=l)
    ax_V.plot(Pos_V ,   T_exp[:-1] ,'k-o',label='Pilote')
    ax_C.plot(   Nl ,   T_exp[ -1] ,'k-o',label='Pilote')
    ax_C.plot([0,Nl],2*[T_exp[ -1]],'k--')

    fig_V.suptitle('Mean Voute Temperature'  ,fontsize=20)
    fig_C.suptitle('Mean Chimney Temperature',fontsize=20)
    fig_V.legend(title='Case',fontsize=16,loc='lower left',bbox_to_anchor=(0.35,0.12))
    ax_V.set_xlabel('X Position [m]'  ,fontsize=16)
    ax_V.set_ylabel('Temperature [K]' ,fontsize=16)
    ax_C.set_xlabel('Cases',fontsize=16)
    ax_C.set_ylabel('Temperature [K]' ,fontsize=16)
    ax_V.set_xticks(Pos_V)
    ax_C.set_xticks(range(Nl+1))
    ax_V.set_xticklabels(['0.51','2.17','3.098'],fontsize=14)
    ax_C.set_xticklabels(labels+['Pilote'],fontsize=14)
    fig_V.subplots_adjust(right=0.5)

    util.SaveFig(fig_V,d_compa+'PLOT/T-Voute.pdf')
    util.SaveFig(fig_C,d_compa+'PLOT/T-Chemine.pdf')