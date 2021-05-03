#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajectory of vortex A from 1 September 2017
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
postdates= True
day1 = datetime(2017,9,1,0)
day2 = day1 + timedelta(days=30)
day1 = datetime(2017,10,1,0)
day2 = datetime(2017,10,18,0)
if predates:
    day1 = datetime(2017,8,24,0)
    day2 = datetime(2017,9,1,0)
    update = False
if postdates:
    day1 = datetime(2017,10,12,0)
    day2 = datetime(2017,10,19,0)
    update = False

if update==False:
# Initial rune
    dats = {}
    i = 0
else:
    # Continuation
    with gzip.open('ERA5-extract-PT5.pkl','rb') as f:
        dats = pickle.load(f)
    i = len(dats)
    print('dats length ',len(dats))

date = day1
pts =   list(np.arange(420,610,5))
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
    if date >= datetime(2017,10,12,0):
        dats[i] = dat.interpolPT(pts,varList=varList,latRange=(5,35),lonRange=(100,220))
    elif date >= datetime(2017,10,3,0):
        datr = dat.shift2west(-180)
        dats[i] = datr.interpolPT(pts,varList=varList,latRange=(20,50),lonRange=(-180,-60))
    elif date >= datetime(2017,9,23,0):
        datr = dat.shift2west(-180)
        dats[i] = datr.interpolPT(pts,varList=varList,latRange=(20,50),lonRange=(-120,0))
    elif  date >= datetime(2017,9,11,0):
        datr = dat.shift2west(-180)
        dats[i] = datr.interpolPT(pts,varList=varList,latRange=(20,50),lonRange=(-30,90))
    elif date >= datetime(2017,9,1,0):
        dats[i] = dat.interpolPT(pts,varList=varList,latRange=(30,60),lonRange=(0,120))
    else:
        datr = dat.shift2west(-180)
        dats[i] = datr.interpolPT(pts,varList=varList,latRange=(40,70),lonRange=(-40,80))
    dat.close()
    date += timedelta(hours=6)
    i += 1

#%%
if predates:
    outbak = 'ERA5-extract-A-PT5-pre.pkl'
    outfile = 'ERA5-extract-A-PT5-new.pkl'
    with gzip.open(outbak,'wb') as f:
        pickle.dump(dats,f)
    with gzip.open('ERA5-extract-A-PT5.pkl','rb') as f:
        dats2 = pickle.load(f)
    for i2 in range(len(dats2)):
        dats[i+i2] = dats2[i2]
elif postdates:
    outfile = 'ERA5-extract-A-PT5-post.pkl'
else:
    outfile = 'ERA5-extract-A-PT5.pkl'
with gzip.open(outfile,'wb') as f:
    pickle.dump(dats,f)
