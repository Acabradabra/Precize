#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#====================> Files
dir0='/mnt/scratch/ZEUS/FLUENT/PRECIZE/RUN-M00/'
# f_param=dir0+'params.tsv'
f_param=dir0+'50_00pct_h2-V00.tsv'

#====================> Combustion parameters
Hyb=0.5 # Hybridation (Hydrogen power fraction)
# Hyb=1 # Hybridation (Hydrogen power fraction)
# Pow=800 # Power (KW)
# Pow=100 # Power (KW)
Pow=10 # Power (KW)
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

#====================> Geometry
d_fc=10.22 # [mm] diameter fuel central
d_fs=11.22 # [mm] diameter fuel side
d_oh=7.0   # [mm] diameter oxygen hublo
d_oc=18.1  # [mm] diameter oxygen central
d_os=10    # [mm] diameter oxygen side
s_bec=14145 # [mm2] section bec
p_bec=475   # [mm] perimeter bec
s_pyro=50   # [mm] size pyro
s_chem=1474875 # [mm2] section chemney
p_chem=8100    # [mm] perimeter chemney
d_vis=154.5 # [mm] diameter visse
h_out=2e3  # [mm] height of the outlet
x_out=1260 # [mm] length of the outlet
y_out=1500 # [mm] length of the outlet

#====================> Klinker parameters
mdot=5 # [T/j] production klinker

y_alooh_feed=0.3400
y_caco3_feed=0.5200
y_humid_feed=0.0075
y_h2o_alooh=0.150156
y_co2_caco3=0.439712