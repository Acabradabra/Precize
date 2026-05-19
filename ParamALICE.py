#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#====================> Files
dfl='/mnt/beegfs/ZEUS/FLUENT/'
dir0='/mnt/beegfs/ZEUS/FLUENT/ALICE/RUN-00/'
dir0=dfl+'ALICE/RUN-00/'
Vs='V00'

#====================> Runs
Dirs=[
'DUMP-01-Ignit/'
]

#====================> Atmospheric conditions
P_atm=101325 # [Pa] Atmospheric pressure
T_atm=300 # [K] Atmospheric temperature
Rho_air=1.17 # [kg/m3] Density of air at atmospheric conditions
# Rho_air=1.15 # [kg/m3] Density of air at atmospheric conditions

#====================> Combustion parameters
Hyb=0 # Hybridation (Hydrogen power fraction)
Imp=0.6 # Impulse (Jet/Anular)
Pow=800 # Power (KW)
if Hyb==1 :
    Rep=0 # Repartition (Top/Bottom)
    l_gn=1    # Excess air ratio Natural gas 
    l_h2=1    # Excess air ratio Hydrogen
    Rb_o=0.15 # Repartition bottom oxygen (central)
elif Hyb==0.5 :
     Rep=0.5 # Repartition (Top/Bottom)
     l_gn=0.95 # Excess air ratio Natural gas 
     l_h2=1    # Excess air ratio Hydrogen
     Rb_o=0.38 # Repartition bottom oxygen (central)
elif Hyb==0 :
     Rep=0.5 # Repartition (Top/Bottom)
     l_gn=1  # Excess air ratio Natural gas 
     l_h2=1    # Excess air ratio Hydrogen
     Rb_o=0.38 # Repartition bottom oxygen (central)

N2  =True # Nitrogen in GN
C3H8=False # Propane in GN

Q_h=0 # [m3/h] Air flow rate for hublo

#====================> Name
f_param=dir0+f'SET/Set_{Hyb*100:.0f}ph2_{Pow*Rep:.0f}kWt_{Pow*(1-Rep):.0f}kWb_{l_gn*100:.0f}lgn_{l_h2*100:.0f}lh2_{Vs}.tsv'

#====================> Geometry
d_fc=10.22 # [mm] diameter fuel central
d_fs=11.22 # [mm] diameter fuel side
d_oh=7.0   # [mm] diameter oxygen hublo
d_oc=18.1  # [mm] diameter oxygen central
d_os=10    # [mm] diameter oxygen side

#====================> Original Parameters
GN_compo={ 'CH4':0.9205,'C2H6':0.0417,'C3H8':0.0093,'C4H10':0.0016,'C5H12':0.0005,'C6H14':0.0032, 'CO2':0.01,'N2':0.0132} # Molar
MF_feeds={'f_hp_l':0.00634,'f_hp_u':0.0038,'f_lp_l':0.00422,'f_lp_u':0.00253,'o_p_l':0.01528,'o_p_u':0.02413,'o_s':2*0.01247} #kg/s

#===============> Param Visualisation
cmesh=0 # coef coordinate (0 : no ticks)
# cmesh=1e2 # coef coordinate (0 : no ticks)
Vars=['Tx']
Param_visu={
'Tx'   :['Data-x0.dat'  ,'yz','temperature'       ,'Temperature [K]'     ,[],[]   ,[],cmesh,[]         ,(5,9) ,'inferno','X0-Temperature.png' ,['ISO',[1800]]],
# 'Tcol' :['Data-F-yc.dat','xz','temperature'       ,'Temperature [K]'     ,[-0.3,0.3],[2,4]   ,[],cmesh,[]         ,(5,9) ,'inferno','Col-Temperature.png' ,['ISO',[1800]]],
# 'Vcol' :['Data-F-yc.dat','xz','velocity-magnitude','Velocity [m/s]'      ,[-0.3,0.3],[2,3]   ,[],cmesh,[]         ,(7,9) ,'cividis','Col-Velocity.png'    ,['QUIV',[2e-2,2e-2,300]]],
# 'YCO2f':['Data-F-x0.dat','yz','co2'               ,'$Y_{CO_2}$ [-]'      ,[]        ,[-0.2,2],[],cmesh,[]         ,(15,7),'viridis','Four-YCO2.png'       ,[]            ],
# 'YH2Of':['Data-F-x0.dat','yz','h2o'               ,'$Y_{H_2O}$ [-]'      ,[]        ,[-0.2,2],[],cmesh,[]         ,(15,7),'viridis','Four-YH2O.png'       ,[]            ],
# 'Tf'   :['Data-F-x0.dat','yz','temperature'       ,'Temperature [K]'     ,[]        ,[-0.2,2],[],cmesh,[ 290,2500],(15,7),'inferno','Four-Temperature.png',['ISO',[2000,2500]]],
# 'Ptf'  :['Data-F-x0.dat','yz','total-pressure'    ,'Total pressure [Pa]' ,[]        ,[-0.2,2],[],cmesh,[]         ,(15,7),'cividis','Four-Pt.png'         ,['ISO',[0]]        ],
# 'Psf'  :['Data-F-x0.dat','yz','pressure'          ,'Static pressure [Pa]',[]        ,[-0.2,2],[],cmesh,[]         ,(15,7),'cividis','Four-Ps.png'         ,['ISO',[0]]        ],
# 'Velf' :['Data-F-x0.dat','yz','velocity-magnitude','Velocity [m/s]'      ,[]        ,[-0.2,2],[],cmesh,[]         ,(15,7),'cividis','Four-Velocity.png'   ,['QUIV',[1e-1,1e-1,300]]],
# 'Ttop' :['Data-TOP.dat' ,'yx','temperature'       ,'Temperature [K]'     ,[]        ,[]      ,[],cmesh,[]         ,(15,7),'inferno','Top-Temperature.png' ,['INTERP']    ],
# 'Htop' :['Data-TOP.dat' ,'yx','heat-flux'         ,'Heat flux [W/m2]'    ,[]        ,[]      ,[],cmesh,[]         ,(15,7),'inferno','Top-HeatFlux.png'    ,['INTERP']    ],
}