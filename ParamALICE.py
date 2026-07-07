#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from numpy import *

#====================> Files
dfl='/mnt/beegfs/ZEUS/FLUENT/'
# dir0=dfl+'ALICE/RUN-00/'
dir0=dfl+'ALICE/RUN-01/'
Vs='V00'

if   dir0[-2]=='0' : pmid='zy' ; pflo='zx'
elif dir0[-2]=='1' : pmid='yz' ; pflo='yx' 

#====================> Runs
Dirs=[
'DUMP-03-DO/'
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
     # Rep=0.5 # Repartition (Top/Bottom)
     Rep=0.375 # Repartition (Top/Bottom)
     l_gn=1.01 # Excess air ratio Natural gas 
     l_h2=1.01 # Excess air ratio Hydrogen
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
# z0=-0.75
z0=-0.575
tol=1e-3

#====================> Boundary conditions
eps_rf=0.80
eps_ss=1
eps_cb=0.9
htc_wb=150
htc_ta=20
htc_tb=50
ctt=0.1
temp_water=300
temp_ext=300
h_wl=0.25 # htc wall left 
h_wr=0.25 # htc wall right
h_wt=0.25 # htc wall top
h_wb=0.30 # htc wall bot 
h_wf=0.35 # htc wall front 
h_cv=0.35 # htc wall convergent
h_cd=0.40 # htc wall condui
h_vo=15   # htc opening volet
h_vw=0.3  # htc wall volet
h_fi=0.3 # htc floor isolated
l_rf=0.48 # [W/m-K] thermal conductivity refractaire

#====================> Original Parameters
GN_compo={ 'CH4':0.9205,'C2H6':0.0417,'C3H8':0.0093,'C4H10':0.0016,'C5H12':0.0005,'C6H14':0.0032, 'CO2':0.01,'N2':0.0132} # Molar
MF_feeds={'f_hp_l':0.00634,'f_hp_u':0.0038,'f_lp_l':0.00422,'f_lp_u':0.00253,'o_p_l':0.01528,'o_p_u':0.02413,'o_s':2*0.01247} #kg/s

#===============> Data Pilote Temperature
Pos_a=arange(0.5,5.6,0.5)
Pos_r=array([-0.6,0,0.6])
Pos_z=array([1,2.5])
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
cmesh=-1 # coef coordinate (-1 : original)
# cmesh=0 # coef coordinate (0 : no ticks)
# cmesh=1 # coef coordinate (1 : rounded)
Vars=['Vx','Vx_z']
Param_visu={
'Vabs' :['DRad-mid.dat'  ,pmid,'volumetric-absorbed-radiation',r'V Iabs [$W/m^3$]',[],[]        , [],cmesh,[]     , (25,5) , 'inferno','Mid-Vabs.png'            ,[]],
'Vemt' :['DRad-mid.dat'  ,pmid,'volumetric-emitted-radiation' ,r'V Iemt [$W/m^3$]',[],[]        , [],cmesh,[]     , (25,5) , 'inferno','Mid-Vemt.png'            ,[]],
'abs'  :['DRad-mid.dat'  ,pmid,'absorption-coefficient',r'Abs [-]'        , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-Abs.png'             ,[]],
'A0'   :['DRad-mid.dat'  ,pmid,'a0'                ,r'Abs gg0 [-]'        , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-A0.png'              ,['Tiso',[300,2400]]],
'A1'   :['DRad-mid.dat'  ,pmid,'a1'                ,r'Abs gg1 [-]'        , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-A1.png'              ,['Tiso',[300,2400]]],
'A2'   :['DRad-mid.dat'  ,pmid,'a2'                ,r'Abs gg2 [-]'        , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-A2.png'              ,['Tiso',[300,2400]]],
'A3'   :['DRad-mid.dat'  ,pmid,'a3'                ,r'Abs gg3 [-]'        , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-A3.png'              ,['Tiso',[300,2400]]],
'A4'   :['DRad-mid.dat'  ,pmid,'a4'                ,r'Abs gg4 [-]'        , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-A4.png'              ,['Tiso',[300,2400]]],
'K0'   :['DRad-mid.dat'  ,pmid,'k0'                ,r'Coef gg0 [-]'       , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-K0.png'              ,['Tiso',[300,2400]]],
'K1'   :['DRad-mid.dat'  ,pmid,'k1'                ,r'Coef gg1 [-]'       , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-K1.png'              ,['Tiso',[300,2400]]],
'K2'   :['DRad-mid.dat'  ,pmid,'k2'                ,r'Coef gg2 [-]'       , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-K2.png'              ,['Tiso',[300,2400]]],
'K3'   :['DRad-mid.dat'  ,pmid,'k3'                ,r'Coef gg3 [-]'       , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-K3.png'              ,['Tiso',[300,2400]]],
'K4'   :['DRad-mid.dat'  ,pmid,'k4'                ,r'Coef gg4 [-]'       , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-K4.png'              ,['Tiso',[300,2400]]],
'Wm'   :['DRad-mid.dat'  ,pmid,'wm'                ,r'Wm [$g/mol$]'       , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-Wm.png'              ,[]],
'Rhcx' :['Data-mid.dat'  ,pmid,'rhc'          ,r'$X_{H_2O}/X_{CO_2}$ [-]' , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-Rhc.png'             ,[]],
'Hr'   :['Data-mid.dat'  ,pmid,'expr:hrr',r'Heat release rate [$mW/m^3$]' , [-0.1,1],[-0.2,0.2] , [],cmesh,[0,1e3], (15,6) , 'inferno','Mid-HeatRelease.png'     ,['GAIN',[1e-6]]],
'Hbot' :['Data-floor.dat',pflo,'heat-flux'         ,r'Heat flux [$W/m^2$]', []      ,[]         , [],cmesh,[]     , (15,5) , 'inferno','Floor-HeatFlux.png'      ,['INTERP','ISO',[1e3] ,'SEL',['z',z0-tol,z0+tol]]],
'Tbot' :['Data-floor.dat',pflo,'temperature'       ,r'Temperature [K]'    , []      ,[]         , [],cmesh,[]     , (15,5) , 'inferno','Floor-Temperature.png'   ,['INTERP','ISO',[2000],'SEL',['z',z0-tol,z0+tol]]],
# 'Htop' :['Data-wall.dat' ,pflo,'heat-flux'         ,r'Heat flux [$W/m^2$]', []      ,[]         , [],cmesh,[]     , (15,5) , 'inferno','Roof-HeatFlux.png'       ,['INTERP','ISO',[1e3] ]],
'Ttop' :['Data-roof.dat' ,pflo,'temperature'       ,r'Temperature [K]'    , []      ,[]         , [],cmesh,[]     , (15,5) , 'inferno','Roof-Temperature.png'    ,['INTERP','ISO',[2000]]],
'Tx_z' :['Data-mid.dat'  ,pmid,'temperature'       ,r'Temperature [K]'    , [-0.1,4],[]         , [],cmesh,[]     , (15,6) , 'inferno','Mid-Temperature_Zoom.png',['ISO',[2000,3000,4000]]],
'Tx'   :['Data-mid.dat'  ,pmid,'temperature'       ,r'Temperature [K]'    , []      ,[]         , [],cmesh,[]     , (25,5) , 'inferno','Mid-Temperature.png'     ,['ISO',[2000,3000,4000]]],
'XCO2x':['Data-mid.dat'  ,pmid,'molef-co2'         ,r'$X_{CO_2}$ [%]'     , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-XCO2.png'            ,['GAIN',[1e2]]          ],
'XH2Ox':['Data-mid.dat'  ,pmid,'molef-h2o'         ,r'$X_{H_2O}$ [%]'     , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-XH2O.png'            ,['GAIN',[1e2]]          ],
'XO2x' :['Data-mid.dat'  ,pmid,'molef-o2'          ,r'$X_{O_2}$ [%]'      , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-XO2.png'             ,['GAIN',[1e2]]          ],
'XN2x' :['Data-mid.dat'  ,pmid,'molef-n2'          ,r'$X_{N_2}$ [%]'      , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-XN2.png'             ,['GAIN',[1e2]]          ],
'YCO2x':['Data-mid.dat'  ,pmid,'co2'               ,r'$Y_{CO_2}$ [%]'     , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-YCO2.png'            ,['GAIN',[1e2]]          ],
'YH2Ox':['Data-mid.dat'  ,pmid,'h2o'               ,r'$Y_{H_2O}$ [%]'     , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-YH2O.png'            ,['GAIN',[1e2]]          ],
'YO2x' :['Data-mid.dat'  ,pmid,'o2'                ,r'$Y_{O_2}$ [%]'      , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-YO2.png'             ,['GAIN',[1e2]]          ],
'YN2x' :['Data-mid.dat'  ,pmid,'n2'                ,r'$Y_{N_2}$ [%]'      , []      ,[]         , [],cmesh,[]     , (25,5) , 'viridis','Mid-YN2.png'             ,['GAIN',[1e2]]          ],
'Psx'  :['Data-mid.dat'  ,pmid,'pressure'          ,r'Pstat rel [Pa]'     , []      ,[]         , [],cmesh,[0,2e2], (25,5) , 'cividis','Mid-Pstat.png'           ,[]],
'Ptx'  :['Data-mid.dat'  ,pmid,'total-pressure'    ,r'Ptot rel [Pa]'      , []      ,[]         , [],cmesh,[0,2e2], (25,5) , 'cividis','Mid-Ptot.png'            ,[]],
'Vy'   :['Data-upper.dat',pflo,'velocity-magnitude',r'Velocity [m/s]'     , []      ,[]         , [],cmesh,[]     , (15,5) , 'cividis','Top-Velocity.png'        ,['QUIV',[10e-2,10e-2,500]]],
'Vy_z' :['Data-upper.dat',pflo,'velocity-magnitude',r'Velocity [m/s]'     , [3,7]   ,[]         , [],cmesh,[0,5]  , (10,6) , 'cividis','Top-Velocity_Zoom.png'   ,['QUIV',[ 5e-2, 5e-2,100]]],
'Vx'   :['Data-mid.dat'  ,pmid,'velocity-magnitude',r'Velocity [m/s]'     , []      ,[]         , [],cmesh,[]     , (25,5) , 'cividis','Mid-Velocity.png'        ,['QUIV',[10e-2,10e-2,750]]],
'Vx_z' :['Data-mid.dat'  ,pmid,'velocity-magnitude',r'Velocity [m/s]'     , [3,7]   ,[]         , [],cmesh,[0,5]  , (15,6) , 'cividis','Mid-Velocity_Zoom.png'   ,['QUIV',[ 5e-2, 5e-2,100]]],
}

#===============> Composition position
probe='outlet'