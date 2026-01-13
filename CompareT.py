#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%=================================================================================
from unittest import case
from IPython import get_ipython

ip = get_ipython()
if ip is not None:
    ip.run_line_magic("load_ext", "autoreload")
    ip.run_line_magic("autoreload", "2")
#%%=================================================================================
#                     Modules
#===================================================================================
import Utilities as util
(Sysa,NSysa,Arg)=util.Parseur([],0,'Arg : ')
(                             [])=Arg

from numpy import *
import os
import sys
import csv
import time
import Fluent as fl

t0=time.time()
(plt,mtp)=util.Plot0()
#%%=================================================================================
#                     Parameters
#===================================================================================

dir0='/mnt/beegfs/ZEUS/FLUENT/PRECIZE/Walter-50pH2-60pJet/'
Dirs=[
dir0+'DUMP-03-A0/' ,
dir0+'DUMP-04-Ray/',
dir0+'DUMP-06-P1/' ,
]

for d in Dirs :
    print('=> '+d)
    dird=d+'DATA/'+'report-meta-temperature-rfile.out'
    Dr=fl.Report_read(dird)
    