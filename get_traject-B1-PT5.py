#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajectory of vortex B1 from 26 August 2017
It select a few variables in a moving window. The data are stored
for all times in a pickle file. Reading this pickle file stores all the data
in memory allowing fast processing of diagnostics such as vortex tracking.

Data are interpolated onto isentropic levels with interval 5K within a selected range.

The data are stored in a dictionary indexed by numbers starting from 0.

Selected variables: 'T','VO','O3','PV','Z','U','V'

predates option allows to add new dates at the beginning of the series, requiring a shift
in the numbering.

@author: Bernard Legras
"""
from datetime import datetime, timedelta
from ECMWF_N import ECMWF
import numpy as np
#from zISA import zISA
import pickle,gzip
#import matplotlib.pyplot as plt

# predates is to add new dates at the begiinning of the file
# it generates a special file which is then merged with
# the existing one in combine_dats.py
predates = False
update = False
day1 = datetime(2017,8,26,0)
day2 = datetime(2017,10,15,0)
if predates:
    day1 = datetime(2020,1,4,6)
    day2 = day1 + timedelta(days=3)
    update = False

if update==False:
# Initial rune
    dats = {}
    i = 0
else:
    # Continuation
    with gzip.open('ERA5-extract-B1-PT5.pkl','rb') as f:
        dats = pickle.load(f)
    i = len(dats)
    print('dats length ',len(dats))

date = day1
pts =   list(np.arange(420,560,5))
varList = ['T','VO','O3','PV','Z','U','V']

while date < day2:
    print(date)
    dat = ECMWF('FULL-EA',date,exp='VOZ')
    dat._get_var('T')
    dat._get_var('VO')
    dat._get_var('O3')
    dat._get_var('U')
    dat._get_var('V')
    dat._mkp()
    dat._mkz()
    dat._mkthet()
    dat._mkpv()
    if date >= datetime(2017,10,9,0):
        datr = dat.shift2west(-180)
        dats[i] = datr.interpolPT(pts,varList=varList,lonRange=(-140,-20),latRange=(30,60))
    elif date >= datetime(2017,9,22,0):
        dats[i] = dat.interpolPT(pts,varList=varList,lonRange=(140,260),latRange=(35,65))
    elif  date >= datetime(2017,9,10,0):
        dats[i] = dat.interpolPT(pts,varList=varList,lonRange=(60,180),latRange=(35,65))
    else:
        datr = dat.shift2west(-180)
        dats[i] = datr.interpolPT(pts,varList=varList,lonRange=(-20,100),latRange=(35,65))
    dat.close()
    date += timedelta(hours=6)
    i += 1

#%%
if predates:
    outfile = 'ERA5-extract-pre.pkl'
else:
    outfile = 'ERA5-extract-B1-PT5.pkl'
with gzip.open(outfile,'wb') as f:
    pickle.dump(dats,f)
