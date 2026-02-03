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
(Sysa,NSysa,Arg)=util.Parseur([],1,'Arg : Params')
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

param=Sysa[0]
if   'Walter' in param : from ParamWalter import *
elif 'Pilot'  in param : from ParamPilot  import *
else                   : raise ValueError('Unknown param file !')

# unit='g/s'
unit='kg/h'

#%%=================================================================================
#                     Process
#===================================================================================
if   unit=='g/s' : fac=1e3
elif unit=='kg/h': fac=3600
else             : fac=1

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
print(f'O2 for H2 : {MO_h*fac:.12f}')
print(f'O2 for GN : {MO_c*fac:.12f}')
print(f'O2 Total  : {MO*fac:.12f}')
print(f'O2 Top    : {MO_t*fac:.12f}')
print(f'O2 Bottom : {MO_b*fac:.12f}')
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
print('\n=> '+util.Col('g','Global flow rates [T0,P0] :'))
print(f'H2 => D0 : {R0_H2:.3f} [Kg/m3]  ,  Q0 : {Q0_H2:.2f} [Nm3/h]')
print(f'GN => D0 : {R0_GN:.3f} [Kg/m3]  ,  Q0 : {Q0_GN:.2f} [Nm3/h]')
print(f'O2 => D0 : {R0_O2:.3f} [Kg/m3]  ,  Q0 : {Q0_O2:.2f} [Nm3/h]')
