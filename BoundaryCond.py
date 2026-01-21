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
(Sysa,NSysa,Arg)=util.Parseur([],0,'Arg : ')
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

Hyb=0.5 # Hybridation (Hydrogen power fraction)
Pow=800 # Power (KW)
Rep=0.5 # Repartition (Top/Bottom)
Imp=0.6 # Impulse (Jet/Anular)

#%%=================================================================================
#                     Process
#===================================================================================

