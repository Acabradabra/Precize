#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#====================> Files
# dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M00/'
dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M01/'
Vs='V03'

#====================> Atmospheric conditions
P_atm=101325 # [Pa] Atmospheric pressure
T_atm=300 # [K] Atmospheric temperature
Rho_air=1.17 # [kg/m3] Density of air at atmospheric conditions
# Rho_air=1.15 # [kg/m3] Density of air at atmospheric conditions

#====================> Darcy
DP_leak=5  # [Pa] Pressure drop through wall
DP_jupe=10 # [Pa] Pressure drop through wall
L_pyro=0.386
L_viss=0.456
L_jupe=0.010
A_viss=18627 # [mm2] Section visse
A_pyro=25    # [mm2] Section pyro
M_leak=0.45*9.5 # [g/s] Target mass flow rate of leaks
M_jupe=50       # [g/s] Target mass flow rate through jupe
nu=1.85e-5    # Pa.s Dynamic viscosity of air at atmospheric conditions (25d Pa)
mu=nu/Rho_air # [m2/s] Kinematic viscosity of air at atmospheric conditions (25d Pa)

#====================> Combustion parameters
Hyb=1 # Hybridation (Hydrogen power fraction)
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

N2  =False # Nitrogen in GN
C3H8=False # Propane in GN

Q_h=2.5 # [m3/h] Air flow rate for hublo

#====================> Name
f_param=dir0+f'SET/Set_{Hyb*100:.0f}ph2_{Pow*Rep:.0f}kWt_{Pow*(1-Rep):.0f}kWb_{l_gn*100:.0f}lgn_{l_h2*100:.0f}lh2_{Vs}.tsv'

#====================> Bath
Tbath=1600 # [°C] 100% H2

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
s_jb=17600 # [mm2] section jupe-bot
p_jb=3520 # [mm] perimeter jupe-bot

#====================> Convection
h_wb=1000 # [W/(m2 K)] Heat transfer coefficient wall-bec
h_nat=10 # [W/(m2 K)] Heat transfer coefficient Natural convection
h_col=20 # [W/(m2 K)] Heat transfer coefficient Natural convection

#====================> Radiation
e_refract=0.85 # Emissivity refractaire
e_shell=0.80 # Emissivity shell
e_ss304=0.3 # Emissivity stainless steel 304L

#====================> Klinker parameters
mdot_kk=6 # [T/j] production klinker
#=====> Walter
# y_alooh_feed=0.3400
# y_caco3_feed=0.5200
# y_humid_feed=0.0075
# y_h2o_alooh=0.150156
# y_co2_caco3=0.439712
#=====> New
y_alooh_feed=0.512
y_caco3_feed=0.488
y_humid_caco3=0.35e-2
y_humid_alooh=0.5e-2
y_h2o_alooh=0.112
y_feo_alooh=0.2323
y_sio_alooh=0.0427
y_tio_alooh=0.0283
y_alo_alooh=0.5668
y_co2_caco3=0.426
y_cao_caco3=0.539
#=====> heat loss
# heat_sink=-187.030729
heat_sink=-292 # [kW]
