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
from numpy import *
from scipy.interpolate import interp1d
(plt,mtp)=util.Plot0()

#%%=================================================================================
#                     Parameters
#===================================================================================

Ri=(0.5*323.9)-8 # [mm] internal radius
Rg=60 # [mm] guillotine hole radius

x0=10
Np=int(1e3)

#%%=================================================================================
#                     Process
#===================================================================================
def y(x,r) : return(sqrt(r**2-x**2))
#===================================================================================
def Ratio(x0,R,Rg,Np) :
    Vx=array(linspace(0,x0,Np))
    Vy=y(Vx,R)

    S0=pi*R**2
    S1=4*trapezoid(Vy,Vx)
    # S0=305**2
    # S1=305*2*x0
    S2=pi*Rg**2
    
    R1=100*S1/S0
    R2=100*S2/S0
    R=R1+R2
    # print(f'=> Opening ratio R1 : {R1:.2f} %')
    # print(f'=> Opening ratio R2 : {R2:.2f} %')
    # print(f'=> Opening ratio R  : {R :.2f} %')
    return(R1,R2,R)
#===================================================================================

Vx=linspace(10,Ri,int(1e2))
# Vx=linspace(10,142.5,int(1e2))
VR=[]
for x1 in Vx :
    R1,R2,R=Ratio(x1,Ri,Rg,Np)
    VR.append(R)
    # print(f'=> x0 : {x1:.3f} mm  ,  R1 : {R1:.2f} %  ,  R2 : {R2:.2f} %  ,  R : {R:.2f} %')
VR=array(VR)
XR=interp1d(VR,Vx)

VR2=[25,50,75,100]
Vx2=XR(VR2)
Vx2[-1]=152.5

fig,ax=plt.subplots(figsize=(10,7))
ax.plot(Vx,VR,'g')
print(f'=> R int : {Ri:.2f} [mm]')
for r,v in zip(VR2,Vx2) :
    print(f'=> R : {r:.2f} %  ,  x0 : {v:.3f} mm')
    ax.plot(2*[v],[25,100],':k')
    ax.text(v-3,102,f'{v:.0f}' ,fontsize=12,color='k')
ax.plot([Vx2[0],Vx2[-2]],[VR2[0],VR2[-2]],'--r')

ax.set_xlabel('x0 [mm]'                        ,fontsize=15)
ax.set_ylabel('R [%]'                          ,fontsize=15)
ax.set_title('Guillotine surface opening ratio',fontsize=15)
util.SaveFig(fig,'Plot/Guillotine.pdf')
