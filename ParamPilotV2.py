#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#====================> Files
# dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M00/'
dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M01/'

#====================> Atmospheric conditions
P_atm=101325 # [Pa] Atmospheric pressure
T_atm=300 # [K] Atmospheric temperature
# Rho_air=1.17 # [kg/m3] Density of air at atmospheric conditions
Rho_air=1.15 # [kg/m3] Density of air at atmospheric conditions

#====================> Darcy
DP=25 # [Pa] Pressure drop through wall
L_pyro=0.386
L_viss=0.456
A_viss=18627 # [mm2] Section visse
A_pyro=25 # [mm2] Section pyro
M_leak=0.45*9.5 # [g/s] Target mass flow rate of leaks
nu=1.4e-4

#====================> Combustion parameters
Hyb=0.5 # Hybridation (Hydrogen power fraction)
# Hyb=1 # Hybridation (Hydrogen power fraction)
# Pow=800 # Power (KW)
# Pow=100 # Power (KW)
# Pow=10 # Power (KW)
Pow=1 # Power (KW)
Rep=0.5 # Repartition (Top/Bottom)
Imp=0.6 # Impulse (Jet/Anular)
l_gn=0.95 # Excess air ratio Natural gas 
# l_gn=1    # Excess air ratio Natural gas 
l_h2=1    # Excess air ratio Hydrogen
# Rb_o=0.38 # Repartition bottom oxygen (central)
Rb_o=0.15 # Repartition bottom oxygen (central)

N2  =False # Nitrogen in GN
C3H8=False # Propane in GN

Q_h=2.5 # [m3/h] Air flow rate for hublo

#====================> Name
f_param=dir0+f'Set_{Hyb*100:.0f}p-h2_{Pow}kW_V00.tsv'

#====================> Geometry
d_fc=10.22 # [mm] diameter fuel central
d_fs=11.22 # [mm] diameter fuel side
d_oh=7.0   # [mm] diameter oxygen hublo
d_oc=18.1  # [mm] diameter oxygen central
d_os=10    # [mm] diameter oxygen side
s_bec=14145 # [mm2] section bec
p_bec=475   # [mm] perimeter bec
s_pyro=50   # [mm] size pyro
d_vis=104   # [mm] diameter visse
d_chem=307.9 # [mm] diameter chemney
s_out=25397 # [mm2] section outlet
p_out=795   # [mm] perimeter outlet
s_f2=101587 # [mm2] section fluid ext inlet
p_f2=2538   # [mm] perimeter fluid ext inlet
h_jup=455   # [mm] height jupe
h_dil=102   # [mm] height dilution
l_dil=430   # [mm] length dilution
h_jeu=17    # [mm] height jeu
l_jeu=335   # [mm] length jeu

#====================> Klinker parameters
mdot=5 # [T/j] production klinker

y_alooh_feed=0.3400
y_caco3_feed=0.5200
y_humid_feed=0.0075
y_h2o_alooh=0.150156
y_co2_caco3=0.439712