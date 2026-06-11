#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from numpy import *

#====================> Files
dfl='/mnt/beegfs/ZEUS/FLUENT/'
dir0='/mnt/beegfs/ZEUS/FLUENT/ALICE/RUN-00/'
dir0=dfl+'ALICE/RUN-00/'
Vs='V00'

#====================> Runs
Dirs=[
'DUMP-05-DO-UDF2/'
]

#====================> Atmospheric conditions
P_atm=101325 # [Pa] Atmospheric pressure
T_atm=300 # [K] Atmospheric temperature
Rho_air=1.17 # [kg/m3] Density of air at atmospheric conditions
# Rho_air=1.15 # [kg/m3] Density of air at atmospheric conditions

#====================> Radiation parameters
Rhc_lim=[0.01,4]
Rhc_tol=1e-10

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

#===============> Data Pilote Temperature
Pos_a=arange(0.5,5.6,0.5)
Pos_r=array([-0.6,0,0.6])
Tax=array([1497,1507,1539,1561,1564,1571,1546,1472,1454,1455,1437])
Tz1=array([1515,1517,1501])
Tz2=array([1576,1571,1551])

#===============> Data Pilote heat flux
Pos_h=array([0.54,0.98,1.42,1.86,2.3,2.74,3.18,3.62,4.06,4.5,4.94,5.38,5.82]) # [m]
Hf=array([201.93,243.67,251.59,281.20,131.45,35.00,81.34,34.34,-4.40,11.24,7.58,2.49,15.91]) # [kW/m2]

#===============> Compo
Compo={
     'O2' :1.05, #[%vol]
     'CO2':38.5, #[%vol]
     'H2O':51  , #[%vol]
     'CO' :0.17, #[mg/MJ]
     'NOx':38  , #[mg/MJ]
}

#===============> Param Visualisation
cmesh=-1 # coef coordinate (0 : no ticks)
# cmesh=0 # coef coordinate (0 : no ticks)
# cmesh=1 # coef coordinate (0 : no ticks)
Vars=[]
# Vars=['Hr']
# Vars=['YCO2x','YH2Ox','YO2x','YN2x']
# Vars=['XCO2x','XH2Ox','XO2x','XN2x']
Param_visu={
'A0'   :['DRad-mid.dat'  ,'zy','a0'                ,r'Abs gg0 [-]'        , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-A0.png'              ,['Tiso',[300,2400]]],
'A1'   :['DRad-mid.dat'  ,'zy','a1'                ,r'Abs gg1 [-]'        , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-A1.png'              ,['Tiso',[300,2400]]],
'A2'   :['DRad-mid.dat'  ,'zy','a2'                ,r'Abs gg2 [-]'        , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-A2.png'              ,['Tiso',[300,2400]]],
'A3'   :['DRad-mid.dat'  ,'zy','a3'                ,r'Abs gg3 [-]'        , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-A3.png'              ,['Tiso',[300,2400]]],
'A4'   :['DRad-mid.dat'  ,'zy','a4'                ,r'Abs gg4 [-]'        , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-A4.png'              ,['Tiso',[300,2400]]],
'K0'   :['DRad-mid.dat'  ,'zy','k0'                ,r'Coef gg0 [-]'       , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-K0.png'              ,['Tiso',[300,2400]]],
'K1'   :['DRad-mid.dat'  ,'zy','k1'                ,r'Coef gg1 [-]'       , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-K1.png'              ,['Tiso',[300,2400]]],
'K2'   :['DRad-mid.dat'  ,'zy','k2'                ,r'Coef gg2 [-]'       , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-K2.png'              ,['Tiso',[300,2400]]],
'K3'   :['DRad-mid.dat'  ,'zy','k3'                ,r'Coef gg3 [-]'       , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-K3.png'              ,['Tiso',[300,2400]]],
'K4'   :['DRad-mid.dat'  ,'zy','k4'                ,r'Coef gg4 [-]'       , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-K4.png'              ,['Tiso',[300,2400]]],
'Wm'   :['DRad-mid.dat'  ,'zy','wm'                ,r'Wm [$g/mol$]'       , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-Wm.png'              ,[]],
'Rhcx' :['Data-mid.dat'  ,'zy','rhc'           ,r'$X_{H_2O}/X_{CO_2}$ [K]', []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-Rhc.png'             ,[]],
'Hr'   :['Data-mid.dat'  ,'zy','expr:hrr',r'Heat release rate [$mW/m^3$]' , [-0.1,1],[-0.2,0.2] , [],cmesh,[0,1e3], (15,6) , 'inferno','Mid-HeatRelease.png'     ,['GAIN',[1e-6]]],
'Ttop' :['Data-roof.dat' ,'zx','temperature'       ,r'Temperature [K]'    , []      ,[]         , [],cmesh,[]     , (15,5) , 'inferno','Roof-Temperature.png'    ,['INTERP','ISO',[2000]]],
'Tx_z' :['Data-mid.dat'  ,'zy','temperature'       ,r'Temperature [K]'    , [-0.1,4],[]         , [],cmesh,[]     , (15,6) , 'inferno','Mid-Temperature_Zoom.png',['ISO',[2000,3000,4000]]],
'Tx'   :['Data-mid.dat'  ,'zy','temperature'       ,r'Temperature [K]'    , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-Temperature.png'     ,['ISO',[2000,3000,4000]]],
'XCO2x':['Data-mid.dat'  ,'zy','molef-co2'         ,r'$X_{CO_2}$ [%]'     , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-XCO2.png'            ,['GAIN',[1e2]]          ],
'XH2Ox':['Data-mid.dat'  ,'zy','molef-h2o'         ,r'$X_{H_2O}$ [%]'     , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-XH2O.png'            ,['GAIN',[1e2]]          ],
'XO2x' :['Data-mid.dat'  ,'zy','molef-o2'          ,r'$X_{O_2}$ [%]'      , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-XO2.png'             ,['GAIN',[1e2]]          ],
'XN2x' :['Data-mid.dat'  ,'zy','molef-n2'          ,r'$X_{N_2}$ [%]'      , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-XN2.png'             ,['GAIN',[1e2]]          ],
'YCO2x':['Data-mid.dat'  ,'zy','co2'               ,r'$Y_{CO_2}$ [%]'     , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-YCO2.png'            ,['GAIN',[1e2]]          ],
'YH2Ox':['Data-mid.dat'  ,'zy','h2o'               ,r'$Y_{H_2O}$ [%]'     , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-YH2O.png'            ,['GAIN',[1e2]]          ],
'YO2x' :['Data-mid.dat'  ,'zy','o2'                ,r'$Y_{O_2}$ [%]'      , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-YO2.png'             ,['GAIN',[1e2]]          ],
'YN2x' :['Data-mid.dat'  ,'zy','n2'                ,r'$Y_{N_2}$ [%]'      , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-YN2.png'             ,['GAIN',[1e2]]          ],
'Vy'   :['Data-upper.dat','zx','velocity-magnitude',r'Velocity [m/s]'     , []      ,[]         , [],cmesh,[]     , (15,5) , 'cividis','Top-Velocity.png'        ,['QUIV',[10e-2,10e-2,500]]],
'Vy_z' :['Data-upper.dat','zx','velocity-magnitude',r'Velocity [m/s]'     , [3,7]   ,[]         , [],cmesh,[0,5]  , (10,6) , 'cividis','Top-Velocity_Zoom.png'   ,['QUIV',[ 5e-2, 5e-2,100]]],
'Vx'   :['Data-mid.dat'  ,'zy','velocity-magnitude',r'Velocity [m/s]'     , []      ,[]         , [],cmesh,[]     , (25,5) , 'cividis','Mid-Velocity.png'        ,['QUIV',[10e-2,10e-2,750]]],
'Vx_z' :['Data-mid.dat'  ,'zy','velocity-magnitude',r'Velocity [m/s]'     , [3,7]   ,[]         , [],cmesh,[0,5]  , (15,6) , 'cividis','Mid-Velocity_Zoom.png'   ,['QUIV',[ 5e-2, 5e-2,100]]],
# 'Ptf'  :['Data-F-x0.dat','yz','total-pressure'    ,'Total pressure [Pa]' ,[]        ,[-0.2,2],[],cmesh,[]         ,(15,7),'cividis','Four-Pt.png'         ,['ISO',[0]]        ],
# 'Psf'  :['Data-F-x0.dat','yz','pressure'          ,'Static pressure [Pa]',[]        ,[-0.2,2],[],cmesh,[]         ,(15,7),'cividis','Four-Ps.png'         ,['ISO',[0]]        ],
# 'Ttop' :['Data-TOP.dat' ,'yx','temperature'       ,'Temperature [K]'     ,[]        ,[]      ,[],cmesh,[]         ,(15,7),'inferno','Top-Temperature.png' ,['INTERP']    ],
# 'Htop' :['Data-TOP.dat' ,'yx','heat-flux'         ,'Heat flux [W/m2]'    ,[]        ,[]      ,[],cmesh,[]         ,(15,7),'inferno','Top-HeatFlux.png'    ,['INTERP']    ],
}
# Vars=list(Param_visu.keys())

#===============> Composition position
probe='outlet'