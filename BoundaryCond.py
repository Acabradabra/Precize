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
(Sysa,NSysa,Arg)=util.Parseur(['File'],0,'Arg : Params')
(                             [ FILE ])=Arg

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

param='Pilot' #Sysa[0]
if   'Walter' in param : from ParamWalter   import *
elif 'Pilot'  in param : from ParamPilotV2  import *
else                   : raise ValueError('Unknown param file !')

# unit='g/s'
unit='kg/hr'

#%%=================================================================================
#                     Combustion
#===================================================================================
if   unit=='g/s'   : fac=1e3
elif unit=='kg/hr' : fac=3600
else               : fac=1

#====================> Fuel composition GN
if C3H8 : Xf_GN=fl.Xf_GN
else    :
    Xf_GN={ k:x for (k,x) in fl.Xf_GN.items() if k!='C3H8' }
    Xf_GN['C2H6']+=fl.Xf_GN['C3H8']
if not N2 : 
    XN2_0=Xf_GN.pop('N2')
    Xf_GN={ k:Xf_GN[k]/(1-XN2_0) for k in Xf_GN.keys() }
Xo_GN=sum([ fl.St_GN[k]*Xf_GN[k] for k in fl.St_GN.keys() if k in Xf_GN.keys() ])
Yf_GN=fl.Conv_XY(Xf_GN)
PCI_GN=sum([ Yf_GN[k]*fl.PCI[k] for k in Yf_GN.keys() if k in fl.PCI.keys() ])

#====================> Mass flow rates
P_h2=   Hyb *Pow*1e3 # Hydrogen power (W)
P_gn=(1-Hyb)*Pow*1e3 # Natural gas power (W)
M_h2=P_h2/fl.PCI['H2'] # Hydrogen mass flow rate (Kg/s)
M_gn=P_gn/PCI_GN       # Natural gas mass flow rate (Kg/s)
M_f=M_h2+M_gn        # Total fuel mass flow rate (Kg/s)
M_t=    Rep *M_f  # Top mass flow rate (Kg/s)
M_b= (1-Rep)*M_f  # Bottom mass flow rate (Kg/s)
M_tc=   Imp *M_t  # Top central jet mass flow rate (Kg/s)
M_ta=(1-Imp)*M_t  # Top annular jet mass flow rate (Kg/s)
M_bc=   Imp *M_b  # Bot central jet mass flow rate (Kg/s)
M_ba=(1-Imp)*M_b  # Bot annular jet mass flow rate (Kg/s)
print('\n=> '+util.Col('r','Fuel mass flow rates [%s] :'%(unit)))
print(f'Top central jet   : {M_tc*fac:.12f}')
print(f'Top annular jet   : {M_ta*fac:.12f}')
print(f'Bot central jet   : {M_bc*fac:.12f}')
print(f'Bot annular jet   : {M_ba*fac:.12f}')

#====================> Fuel composition
Y_h2=M_h2/(M_h2+M_gn) # Hydrogen mass fraction
Y_gn=1-Y_h2           # Natural gas mass fraction
Y_f={ k:Y_gn*Yf_GN[k] for k in Yf_GN.keys() }
Y_f['H2']=Y_h2
X_f=fl.Conv_YX(Y_f)
print(f'\n=> '+util.Col('r','Fuel composition X : '),sum([ X_f[k] for k in X_f.keys() ]))
for k in X_f.keys() : print(f'=> X {k} : {X_f[k]:.12e}')

#====================> Oxygen flow rate
W_GN=fl.W_moy(Xf_GN)         # Natural gas molar mass (g/mol)
sh=8                         # Stoichiometric coefficient Hydrogen
sc=fl.Mol_m['O2']*Xo_GN/W_GN # Stoichiometric coefficient Natural gas
MO_h=sh*l_h2*M_h2 # Hydrogen oxygen mass flow rate (Kg/s)
MO_c=sc*l_gn*M_gn # Natural gas oxygen mass flow rate (Kg/s)
MO=MO_h+MO_c      # Total oxygen mass flow rate (Kg/s)
MO_t=    Rep *MO  # Top oxygen mass flow rate (Kg/s)
MO_b= (1-Rep)*MO  # Bottom oxygen mass flow rate (Kg/s)
MO_bc=    Rb_o*MO_b # Bottom center oxygen mass flow rate (Kg/s)
MO_bs=(1-Rb_o)*MO_b # Bottom side oxygen mass flow rate (Kg/s)
print('\n=> '+util.Col('b','Oxygen mass flow rates [%s] :'%(unit)))
# print(f'O2 for H2 : {MO_h*fac:.12f}')
# print(f'O2 for GN : {MO_c*fac:.12f}')
# print(f'O2 Total  : {MO*fac:.12f}')
print(f'O2 Top    : {MO_t*fac:.12f}')
# print(f'O2 Bottom : {MO_b*fac:.12f}')
print(f'O2 Central : {MO_bc*fac:.12f}')
print(f'O2 Side    : {MO_bs*fac:.12f}')

#====================> Summary
PRT0=fl.P0/(fl.R*fl.T0)
R0_H2=fl.Mol_m['H2']*PRT0*1e-3 # density [Kg/m3]
R0_GN=          W_GN*PRT0*1e-3 # density [Kg/m3]
R0_O2=fl.Mol_m['O2']*PRT0*1e-3 # density [Kg/m3]
Q0_H2=(3.6e3*M_h2)/R0_H2
Q0_GN=(3.6e3*M_gn)/R0_GN
Q0_O2=(3.6e3*MO  )/R0_O2
MO_hb=Q_h*R0_O2
print('\n=> '+util.Col('g','Global flow rates [T0,P0] :'))
print(f'H2 => D0 : {R0_H2:.3f} [Kg/m3]  ,  Q0 : {Q0_H2:.2f} [Nm3/h]')
print(f'GN => D0 : {R0_GN:.3f} [Kg/m3]  ,  Q0 : {Q0_GN:.2f} [Nm3/h]')
print(f'O2 => D0 : {R0_O2:.3f} [Kg/m3]  ,  Q0 : {Q0_O2:.2f} [Nm3/h]')
print(f'Hublo => Q0 : {Q_h:.2f} [Nm3/h]  ,  MO : {MO_hb:.2f} [Kg/h]')

#%%=================================================================================
#                     Geometry
#===================================================================================

s_fc=pi*(d_fc/2)**2      # Central fuel jet section (mm2)
p_fc=pi*d_fc             # Central fuel jet perimeter (mm)
hd_fc=4*s_fc/p_fc        # Central fuel jet hydraulic diameter (mm)
s_fs=pi*(d_fs/2)**2-s_fc # Side fuel jet section (mm2)
p_fs=pi*d_fs       +p_fc # Side fuel jet perimeter (mm)
hd_fs=4*s_fs/p_fs        # Side fuel jet hydraulic diameter (mm)
s_oh=pi*(d_oh/2)**2      # Oxygen hublo section (mm2)
p_oh=pi*d_oh             # Oxygen hublo perimeter (mm)
hd_oh=4*s_oh/p_oh        # Oxygen hublo hydraulic diameter (mm)
s_oc=pi*(d_oc/2)**2-pi*(d_fs/2)**2 # Oxygen central jet section (mm2)
p_oc=pi*d_oc       +pi*d_fs        # Oxygen central jet perimeter (mm)
hd_oc=4*s_oc/p_oc        # Oxygen central jet hydraulic diameter (mm)
s_os=pi*(d_os/2)**2      # Oxygen side jet section (mm2)
p_os=pi*d_os             # Oxygen side jet perimeter (mm)
hd_os=4*s_os/p_os        # Oxygen side jet hydraulic diameter (mm)

hd_bec=4*s_bec/p_bec      # Hydraulic diameter bec (mm)
hd_pyro=s_pyro            # Hydraulic diameter pyro (mm)
hd_vis =d_vis             # Hydraulic diameter visse (mm)
hd_chem=d_chem            # Hydraulic diameter chemney (mm)
hd_out=4*s_out/p_out      # Hydraulic diameter fluid ext outlet (mm)
hd_f2=4*s_f2/p_f2         # Hydraulic diameter fluid ext inlet (mm)
hd_jb=4*s_jb/p_jb         # Hydraulic diameter fluid ext inlet (mm)
hd_dil=2*h_dil*l_dil/(h_dil+l_dil) # Hydraulic diameter dilution (mm)
hd_jup=2*h_jup*l_dil/(h_jup+l_dil) # Hydraulic diameter jupe (mm)
hd_jeu=2*h_jeu*l_jeu/(h_jeu+l_jeu) # Hydraulic diameter jeu (mm)
# hd_chem=4*s_chem/p_chem   # Hydraulic diameter chemney (mm)
# hd_out_x=2*y_out*h_out/(y_out+h_out) # Hydraulic diameter outlet side (mm)
# hd_out_y=2*x_out*h_out/(x_out+h_out) # Hydraulic diameter outlet side (mm)

#%%=================================================================================
#                     Klinker
#===================================================================================
y_vap_feed=(y_humid_caco3*y_caco3_feed+y_humid_alooh*y_alooh_feed) # mass fraction of humidity in raw material
y_h2o_feed=(                             y_h2o_alooh*y_alooh_feed) # mass fraction of humidity in raw material
y_co2_feed=(y_co2_caco3*y_caco3_feed) # mass fraction of CO2      in raw material
y_kk=1-(y_h2o_feed+y_co2_feed) # mass fraction of klinker in raw material
mdot_rm=mdot_kk/y_kk           # mass flow rate of raw material

print('=> y_vap_feed : ',y_vap_feed)
print('=> y_h2o_feed : ',y_h2o_feed)
print('=> y_co2_feed : ',y_co2_feed)
mdot_h2o=(y_h2o_feed+y_vap_feed)*mdot_rm
mdot_co2= y_co2_feed            *mdot_rm

mdot_h2o*=(1e3/(3600*24)) # [kg/s]
mdot_co2*=(1e3/(3600*24)) # [kg/s]

print('\n=> '+util.Col('r','Mass flux from solid [Kg/s] :'))
print(f'=> H2O : {mdot_h2o*fac:.3f} [{unit:s}]  ,  {mdot_h2o*1e3:.3f} [g/s]')
print(f'=> CO2 : {mdot_co2*fac:.3f} [{unit:s}]  ,  {mdot_co2*1e3:.3f} [g/s]')
print(f'=> Ratio : {mdot_co2/mdot_h2o:.2f} [-]')
print(f'=> KK/RM : {100*y_kk:.2f} [%]')

#%%=================================================================================
#                     Pressure Loss
#===================================================================================
A_tot =(A_viss+3*A_pyro)*1e-6
A_jupe=(4*h_jup*l_dil  )*1e-6
m_leak=1e-3*M_leak/A_tot
m_jupe=1e-3*M_jupe/A_jupe
K_viss=mu*L_viss*m_leak/DP_leak
K_pyro=mu*L_pyro*m_leak/DP_leak
K_jupe=mu*L_jupe*m_jupe/DP_jupe
Rv_viss=1/K_viss
Rv_pyro=1/K_pyro
Rv_jupe=1/K_jupe

print('\n=> '+util.Col('r','Permeability :'))
print(f'=> Permeability viss : {K_viss:.3e} [m2]  =>  Visc res : {Rv_viss:.3e} [1/m2] ')
print(f'=> Permeability pyro : {K_pyro:.3e} [m2]  =>  Visc res : {Rv_pyro:.3e} [1/m2] ')
print(f'=> Permeability jupe : {K_jupe:.3e} [m2]  =>  Visc res : {Rv_jupe:.3e} [1/m2] ')

#%%=================================================================================
#                     Writing
#===================================================================================
if FILE :
    print('\n'+util.Col('b','=> Writing : '+f_param))
    Param={}
    Param["hrr"                   ]=[ 'HeatofReaction/CellVolume' , ''       , "hrr"                   , ''                           ]
    Param["pstat_ext"             ]=[ '-rho_air*g*Position.z'     , ''       , "pstat_ext"             , 'pressure'                   ]
    Param["pdyn_ext"              ]=[ 'pstat_ext+DynamicPressure' , ''       , "pdyn_ext"              , 'pressure'                   ]
    Param["area_tc"               ]=[ "Area(['f-tc'])"            , ''       , "area_tc"               , ''                           ]
    Param["eps_refract"           ]=[ f'{e_refract   :.2f}'  , ''             , "eps_refract"           , ''                          ]
    Param["eps_shell"             ]=[ f'{e_shell     :.2f}'  , ''             , "eps_shell"             , ''                          ]
    Param["eps_ss304"             ]=[ f'{e_ss304     :.2f}'  , ''             , "eps_ss304"             , ''                          ]
    Param["htc_wb"                ]=[ f'{h_wb        :.0f}'  , ' [W/(m^2 K)]' , "htc_wb"                , 'heat-transfer-coefficient' ]
    Param["htc_nat"               ]=[ f'{h_nat       :.0f}'  , ' [W/(m^2 K)]' , "htc_nat"               , 'heat-transfer-coefficient' ]
    Param["htc_col"               ]=[ f'{h_col       :.0f}'  , ' [W/(m^2 K)]' , "htc_col"               , 'heat-transfer-coefficient' ]
    Param["rho_air"               ]=[ f'{Rho_air     :.2f}'  , ' [kg/m^3]'    , "rho_air"               , 'density'                   ]
    Param["mdot_fuel_top_annular" ]=[ f'{M_ta *fac   :.12e}' , ' [%s]'%(unit) , "mdot_fuel_top_annular" , 'mass-flow'                 ]
    Param["mdot_fuel_top_axial"   ]=[ f'{M_tc *fac   :.12e}' , ' [%s]'%(unit) , "mdot_fuel_top_axial"   , 'mass-flow'                 ]
    Param["mdot_fuel_bot_annular" ]=[ f'{M_ba *fac   :.12e}' , ' [%s]'%(unit) , "mdot_fuel_bot_annular" , 'mass-flow'                 ]
    Param["mdot_fuel_bot_axial"   ]=[ f'{M_bc *fac   :.12e}' , ' [%s]'%(unit) , "mdot_fuel_bot_axial"   , 'mass-flow'                 ]
    Param["mdot_oxid_porthole"    ]=[ f'{MO_hb       :.12e}' , ' [%s]'%(unit) , "mdot_oxid_porthole"    , 'mass-flow'                 ]
    Param["mdot_oxid_top"         ]=[ f'{MO_t *fac   :.12e}' , ' [%s]'%(unit) , "mdot_oxid_top"         , 'mass-flow'                 ]
    Param["mdot_oxid_bottom"      ]=[ f'{MO_bc*fac   :.12e}' , ' [%s]'%(unit) , "mdot_oxid_bottom"      , 'mass-flow'                 ]
    Param["mdot_oxid_staged"      ]=[ f'{MO_bs*fac   :.12e}' , ' [%s]'%(unit) , "mdot_oxid_staged"      , 'mass-flow'                 ]
    Param["mdot_h2o"              ]=[ f'{mdot_h2o*fac:.12e}' , ' [%s]'%(unit) , "mdot_h2o"              , 'mass-flow'                 ]
    Param["mdot_co2"              ]=[ f'{mdot_co2*fac:.12e}' , ' [%s]'%(unit) , "mdot_co2"              , 'mass-flow'                 ]
    Param["heat_sink"             ]=[ f'{heat_sink   :.6f}'  , ' [kW]'        , "heat_sink"             , ''                          ]
    Param["smdot_h2o"             ]=[ 'mdot_h2o/area_tc'     , ''             , "smdot_h2o"             , ''                          ]
    Param["smdot_co2"             ]=[ 'mdot_co2/area_tc'     , ''             , "smdot_co2"             , ''                          ]
    Param["sheat_sink"            ]=[ 'heat_sink/area_tc'    , ''             , "sheat_sink"            , ''                          ]
    Param["temp_bath"             ]=[ f'{Tbath+273.15:.0f}'  , ' [K]'         , "temp_bath"             , 'temperature'               ]
    Param["temp_ext"              ]=[ f'{300         :.0f}'  , ' [K]'         , "temp_ext"              , 'temperature'               ]
    Param["temp_fuel"             ]=[ f'{300         :.0f}'  , ' [K]'         , "temp_fuel"             , 'temperature'               ]
    Param["temp_oxid"             ]=[ f'{300         :.0f}'  , ' [K]'         , "temp_oxid"             , 'temperature'               ]
    Param["temp_water"            ]=[ f'{300         :.0f}'  , ' [K]'         , "temp_water"            , 'temperature'               ]
    Param["x_fuel_c2h6"           ]=[ f'{X_f['C2H6'] :.12e}' , ''             , "x_fuel_c2h6"           , ''                          ]
    Param["x_fuel_ch4"            ]=[ f'{X_f['CH4' ] :.12e}' , ''             , "x_fuel_ch4"            , ''                          ]
    Param["x_fuel_co2"            ]=[ f'{X_f['CO2' ] :.12e}' , ''             , "x_fuel_co2"            , ''                          ]
    Param["x_fuel_h2"             ]=[ f'{X_f['H2'  ] :.12e}' , ''             , "x_fuel_h2"             , ''                          ]
    Param["x_oxid_o2"             ]=[ f'{ 1          :.0f}'  , ''             , "x_oxid_o2"             , ''                          ]
    Param["hd_fc"                 ]=[ f'{ hd_fc      :.3f}'  , ' [mm]'        , "hd_fc"                 , 'length'                    ]
    Param["hd_fs"                 ]=[ f'{ hd_fs      :.3f}'  , ' [mm]'        , "hd_fs"                 , 'length'                    ]
    Param["hd_oh"                 ]=[ f'{ hd_oh      :.3f}'  , ' [mm]'        , "hd_oh"                 , 'length'                    ]
    Param["hd_oc"                 ]=[ f'{ hd_oc      :.3f}'  , ' [mm]'        , "hd_oc"                 , 'length'                    ]
    Param["hd_os"                 ]=[ f'{ hd_os      :.3f}'  , ' [mm]'        , "hd_os"                 , 'length'                    ]
    Param["hd_bec"                ]=[ f'{ hd_bec     :.3f}'  , ' [mm]'        , "hd_bec"                , 'length'                    ]
    Param["hd_pyro"               ]=[ f'{ hd_pyro    :.3f}'  , ' [mm]'        , "hd_pyro"               , 'length'                    ]
    Param["hd_vis"                ]=[ f'{ hd_vis     :.3f}'  , ' [mm]'        , "hd_vis"                , 'length'                    ]
    Param["hd_chem"               ]=[ f'{ hd_chem    :.3f}'  , ' [mm]'        , "hd_chem"               , 'length'                    ]
    Param["hd_out"                ]=[ f'{ hd_out     :.3f}'  , ' [mm]'        , "hd_out"                , 'length'                    ]
    Param["hd_f2"                 ]=[ f'{ hd_f2      :.3f}'  , ' [mm]'        , "hd_f2"                 , 'length'                    ]
    Param["hd_jb"                 ]=[ f'{ hd_jb      :.3f}'  , ' [mm]'        , "hd_jb"                 , 'length'                    ]
    Param["hd_dil"                ]=[ f'{ hd_dil     :.3f}'  , ' [mm]'        , "hd_dil"                , 'length'                    ]
    Param["hd_jup"                ]=[ f'{ hd_jup     :.3f}'  , ' [mm]'        , "hd_jup"                , 'length'                    ]
    Param["hd_jeu"                ]=[ f'{ hd_jeu     :.3f}'  , ' [mm]'        , "hd_jeu"                , 'length'                    ]
    Param["rv_viss"               ]=[ f'{ Rv_viss    :.3e}'  , ' [m^-2]'      , "rv_viss"               , ''                          ]
    Param["rv_pyro"               ]=[ f'{ Rv_pyro    :.3e}'  , ' [m^-2]'      , "rv_pyro"               , ''                          ]
    Param["rv_jupe"               ]=[ f'{ Rv_jupe    :.3e}'  , ' [m^-2]'      , "rv_jupe"               , ''                          ]

    with open(f_param,'w') as f :
        f.write( 'name\tdefinition\tdescription\tparameterid\tparametername\tunit\tinput-parameter\toutput-parameter\t\n')
        for k in Param.keys() : f.write(f'"{k}"\t"{Param[k][0]+Param[k][1]}"\t""\t""\t"{Param[k][2]}"\t"{Param[k][3]}"\t#f\t#f\t\n')
    f.closed