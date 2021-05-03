#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajectory of Canadian vortex B1
Localize the vortex
Print the positions
Calculate the total path.

Created on Friday 25 September 2020

@author: Bernard Legras
"""
from datetime import datetime, timedelta
#from ECMWF_N import ECMWF
import numpy as np
#from zISA import zISA
import constants as cst
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D # yes it is used
#from matplotlib import cm
#from matplotlib.text import TextPath
import gzip,pickle
import socket
import cartopy.crs as ccrs
#import deepdish as dd
from os.path import join
#from PIL import Image, ImageDraw, ImageFont
#from os import chdir
#from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
import imageio

#%%
pts =   list(np.arange(420,560,5))
LaitFactor = (np.array(pts)/500)**(-4.5)

def tracker(dats,var,lon,lat,lower=0,upper=19,idx=6,jdy=3):
    # extract volume surrounding the target point
    jy = np.where(dats.attr['lats']>=lat)[0][0]
    ix = np.where(dats.attr['lons']>=lon)[0][0]
    #print('jy, ix',jy,ix)
    jmin = max(jy-jdy,0)
    imin = max(ix-idx,0)
    jmax = min(jy+jdy+1,dats.nlat)
    imax = min(ix+idx+1,dats.nlon)
    sample = dats.var[var][lower:upper,jmin:jmax,imin:imax]
    # find the 3d index of the max vorticity in the sample cube
    aa = np.unravel_index(np.argmin(sample, axis=None), sample.shape)
    # return [lat,lon,theta,p,z,PV,vo,T,O3,kz,jz,iz]
    return([dats.attr['lats'][jmin+aa[1]],
            dats.attr['lons'][imin+aa[2]],
            pts[lower+aa[0]],
            cst.p0 * (dats.var['T'][lower+aa[0],jmin+aa[1],imin+aa[2]]/pts[lower+aa[0]])**(1/cst.kappa),
            dats.var['Z'][lower+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['PV'][lower+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['PV'][lower+aa[0],jmin+aa[1],imin+aa[2]]-np.mean(dats.var['PV'][lower+aa[0],jmin+aa[1],:]),
            dats.var['VO'][lower+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['T'][lower+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['O3'][lower+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['O3prime'][lower+aa[0],jmin+aa[1],imin+aa[2]],
            lower+aa[0],jmin+aa[1],imin+aa[2]])

if 'gort' == socket.gethostname():
    rootdir = '/dkol/data/STC/STC-CAN'
elif 'satie' in socket.gethostname():
    rootdir = '/data/STC/STC-CAN'

figsav = False
figargs = dict(bbox_inches='tight',dpi=300)
#%%
with gzip.open('ERA5-extract-B1-PT5.pkl','rb') as f:
    dats = pickle.load(f)
print(len(dats))

#%%
figsav = False
figargs = dict(bbox_inches='tight',dpi=300)
# tracking of the 3D position of the vortex
trac={'dates':[],'lons':[],'lats':[],'vo':[],'z':[],'T':[],'p':[],
      'pv':[],'pvano':[],'pt':[],'o3':[],'o3prime':[],'ix':[],'jy':[],'kz':[]}
# initial position
#lon = 245
#lat = -54
#date = datetime(2020,1,7,6)
lat = 48
lon = -5
upper = 10
lower = 8
date=datetime(2017,8,27,12)
ioff = 6
for i in range(len(dats)-ioff):
    print(i)
    i1 = i + ioff
    idx = 6
    jdy = 2
    var = 'O3prime'
    var = 'LPV'
    # track with one iteration
    try:
        lon = 2*lon - trac['lons'][-2]
        lat = 2*lat - trac['lats'][-2]
        upper = min(28,trac['kz'][-1] + 3)
        lower = max(8,trac['kz'][-1] - 2)
    except: pass
    if (var == 'LPV'):
        if i == 164: lon = 232
        if i == 165: lon = 235
        if i == 166:
            lon = 237; lat = 44
        if i == 167:
            lon = 239; lat = 43
        if i == 168:
            lon = 239; lat = 42
        if i == 169:
            lon = 242; lat = 41
    if (var == 'O3prime'):
        pass

    # change of longitude origin
    if i == 170: lon -= 360
    if i == 171: lon += 360

    dats[i1].var['LPV'] = dats[i1].var['PV']*LaitFactor[:,None,None]
    meanO3 = np.mean(dats[i].var['O3'],axis=(1,2))
    #dats[i1].var['O3prime'] = dats[i1].var['O3']/meanO3[:,None,None] -1
    dats[i1].var['O3prime'] = dats[i1].var['O3']- meanO3[:,None,None]

    try:
        [lat,lon,pt,p,z,pv,pvano,vo,T,o3,o3prime,kz,jy,ix] = \
            tracker(dats[i1],var,lon,lat,upper=upper,lower=lower,idx=idx,jdy=jdy)
        #[lat,lon,pt,p,z,pv,vo,T,o3,kz,jy,ix] = tracker(dats[i],var,lon,lat,upper=upper,lower=lower,idx=idx,jdy=jdy)
    except:
        print('tracking error at step ',i)
        print('tracking terminated')
        break
    trac['dates'].append(date)
    trac['lons'].append(lon)
    trac['lats'].append(lat)
    trac['vo'].append(vo)
    trac['pv'].append(pv)
    trac['pvano'].append(pvano)
    trac['z'].append(z/1000)
    trac['T'].append(T)
    trac['o3'].append(o3)
    trac['o3prime'].append(o3prime)
    trac['p'].append(p)
    trac['pt'].append(T*(cst.p0/p)**cst.kappa)
    trac['kz'].append(kz)
    trac['jy'].append(jy)
    trac['ix'].append(ix)
    date += timedelta(hours=6)

if len(trac['dates']) == len(dats)-ioff:
    print('CONGRATULATIONS, YOU TRACKED THE VORTEX UNTIL THE END')
else:
    print('Tracking performed until ',trac['dates'][-1])
if var=='LPV':
    trac_PV = {}
    for k in trac.keys():
        trac_PV[k] = trac[k][:194].copy()
    trac1 = trac_PV
elif var == 'O3prime':
    trac_O3 = {}
    for k in trac.keys():
        trac_O3[k] = trac[k][:194].copy()
    trac2 = trac_O3

#%% print the positions as a function of time
for i in range(len(trac1['dates'])):
    # kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    print(i,trac1['dates'][i],trac1['lons'][i],trac1['lats'][i],'{:2.1f}'.format(trac1['z'][i]),trac1['kz'][i])
    print(i,trac2['dates'][i],trac2['lons'][i],trac2['lats'][i],'{:2.1f}'.format(trac2['z'][i]),trac2['kz'][i])
# beware that saving here will loose the wind tracking made in wind-census
pickle.dump([trac_PV,trac_O3],open('VortexB1-track.pkl','wb'))

#%% Total displacement
[trac1,trac2] = pickle.load(open('VortexB1-track.pkl','rb'))
# with PV
lats = np.array(trac1['lats'])
lons = np.array(trac1['lons'])
lons[lons<-50] += 360
dy = lats[1:]-lats[:-1]
dx = (lons[1:]-lons[:-1])*np.cos(np.deg2rad(0.5*(lats[1:]+lats[:-1])))
ds = np.sqrt(dx**2+dy**2)
print('total path PV ',(2*np.pi*6371/360)*np.sum(ds))
# with O3
lats = np.array(trac2['lats'])
lons = np.array(trac2['lons'])
lons[lons<-50] += 360
dy = lats[1:]-lats[:-1]
dx = (lons[1:]-lons[:-1])*np.cos(np.deg2rad(0.5*(lats[1:]+lats[:-1])))
ds = np.sqrt(dx**2+dy**2)
print('total path O3 ',(2*np.pi*6371/360)*np.sum(ds))

#%% Checker

#for i in range(50,len(trac1['dates'])):
images = []
movieOn = True
for i in range(0,194):
    # generates scaled Montogomery potential
    Mo = (cst.g * dats[i+ioff].var['Z'] + cst.Cp * dats[i+ioff].var['T']) * 1.e-3
    kz1 = trac1['kz'][i]
    kz2 = trac2['kz'][i]
    fig = plt.figure(figsize=(11,8))
    if dats[i+ioff].attr['lons'][-1] > 180: cm_lon=180
    else: cm_lon=0
    projplate = ccrs.PlateCarree(central_longitude=cm_lon)
    gs = fig.add_gridspec(2,1)
    ax1 = fig.add_subplot(gs[0,0],projection=projplate)
    ax1 = dats[i+ioff].show('PV',kz1,figsize=None,axf=ax1,show=False,scale=1.e6,xylim=True,
                 txt='PV (PVU)  '+str(i)+'  '+str(kz1)+trac1['dates'][i].strftime('  %Y %b %d %HUTC  ')+str(pts[kz1])+'K'+
                 '   lon'+str(trac1['lons'][i])+'   lat'+str(trac1['lats'][i]))
    ax1.plot(trac1['lons'][i]-cm_lon,trac1['lats'][i],'black',marker='x',ms=21,mew=4)
    Mof = gaussian_filter(Mo[kz1,...],4)
    Mof = Mof - np.mean(Mof)
    cs = ax1.contour(Mof,transform=projplate,extent=(dats[i+ioff].attr['lons'][0]-cm_lon,
                dats[i+ioff].attr['lons'][-1]-cm_lon,dats[i+ioff].attr['lats'][0],dats[i+ioff].attr['lats'][-1]),origin='lower')
    plt.clabel(cs)
    ax2 = fig.add_subplot(gs[1,0],projection=projplate)
    ax2 = dats[i+ioff].show('O3',kz2,figsize=None,axf=ax2,show=False,scale=1.e6,xylim=True,
                 txt='O3 (ppmv)  '+str(i)+'  '+str(kz2)+trac2['dates'][i].strftime('  %Y %b %d %HUTC  ')+str(pts[kz2])+'K'+
                 '   lon'+str(trac2['lons'][i])+'   lat'+str(trac2['lats'][i]))
    ax2.plot(trac2['lons'][i]-cm_lon,trac2['lats'][i],'black',marker='x',ms=21,mew=4)
    Mof = gaussian_filter(Mo[kz2,...],4)
    Mof = Mof - np.mean(Mof)
    cs = ax2.contour(Mof,transform=projplate,extent=(dats[i+ioff].attr['lons'][0]-cm_lon,
                dats[i+ioff].attr['lons'][-1]-cm_lon,dats[i+ioff].attr['lats'][0],dats[i+ioff].attr['lats'][-1]),origin='lower')
    plt.clabel(cs)
    ax1.text(-0.05,1.05,'a)',transform=ax1.transAxes,fontsize=22)
    ax2.text(-0.05,1.05,'b)',transform=ax2.transAxes,fontsize=22)
    if movieOn:
        plt.savefig('tmpimB1.png',bbox_inches='tight')
        images.append(imageio.imread('tmpimB1.png'))
    plt.show()
#%%
if movieOn:
    imageio.mimsave('movie_VortexB1.gif', images, fps=2)
    # reshape images to have the same size (required by mp4), the problem occurs in lon
    # find the min size
    minX = 20000
    for j in range(len(images)):
        minX = min(minX,images[j].shape[1])
    for j in range(len(images)):
        if images[j].shape[1] > minX:
            images[j] = images[j][:,:minX,:]

    imageio.mimsave('movie_VortexB1.mp4', images, fps=2)

#%% Total displacement
# trac['lats'] = np.array(trac['lats'])
# trac['lons'] = np.array(trac['lons'])
# dy = trac['lats'][1:]-trac['lats'][:-1]
# dx = (trac['lons'][1:]-trac['lons'][:-1])*np.cos(np.deg2rad(0.5*(trac['lats'][1:]+trac['lats'][:-1])))
# # Correction for crossing Greenwich
# dx[149] = ((trac['lons'][150]-trac['lons'][149]%360))*np.cos(np.deg2rad(0.5*(trac['lats'][149]+trac['lats'][150])))
# ds = np.sqrt(dx**2+dy**2)
# print('total path ',(2*np.pi*6371/360)*np.sum(ds))

#%% Localisation of the vortex

fig,((ax0,ax1),(ax2,ax3),(ax4,ax5))=plt.subplots(3,2,figsize=(10,10),sharex=True)
pvl1 = np.array(trac1['pvano'])*(np.array(trac1['pt'])/500)**(-4.5)
pvl2 = np.array(trac2['pvano'])*(np.array(trac2['pt'])/500)**(-4.5)
ax0.plot(
         trac2['dates'],trac2['lats'],
         trac1['dates'],trac1['lats'],linewidth=2)
trac1l = np.array(trac1['lons'])
trac2l = np.array(trac2['lons'])
trac1l[170:] += 360
trac2l[170:] += 360
ax1.plot(
         trac2['dates'],trac2l,
         trac1['dates'],trac1l,linewidth=2)
ax2.plot(
         trac2['dates'],trac2['pt'],
         trac1['dates'],trac1['pt'],
         linewidth=2)
ax2.legend(loc='upper left',labels=('O3 anomaly','Lait PV'))
# ax3.plot(trac3['dates'],1.e6*np.array(trac3['pvano']),
#          trac2['dates'],1.e6*np.array(trac2['pvano']),
#          trac1['dates'],1.e6*np.array(trac1['pvano']),linewidth=2)
ax3.plot(
         trac2['dates'],1.e6*pvl2,
         trac1['dates'],1.e6*pvl1,linewidth=2)
ax4.plot(
         trac2['dates'],1.e6*np.array(trac2['o3']),
         trac1['dates'],1.e6*np.array(trac1['o3']),linewidth=2)
ax5.plot(
         trac2['dates'],trac2['o3prime'],
         trac1['dates'],trac1['o3prime'],linewidth=2)
ax0.set_ylabel('Latitude')
ax1.set_ylabel('Longitude')
ax2.set_ylabel('PT altitude (K)')
ax3.set_ylabel('Lait PV anomaly(PVU)')
ax4.set_ylabel('Ozone (ppmm)')
ax5.set_ylabel('Normalized ozone anomaly')
fig.suptitle('Canadian event   vortex B1 tracking')
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('VortexMotion_VortexB1_ERA5PV.png'),**figargs)
plt.show()

#%% Plot the horizontal displacement
# plt.plot(trac['lons'],trac['lats'],linewidth=3)
# plt.xlabel('Longitude')
# plt.ylabel('Latitude')
# plt.title('Horizontal displacement between 07-01-2020 and 20-02-2020')
# #mm = {6:'7 Jan',14:'11 Jan',22:'15 Jan',32:'20 Jan',42:'25 Jan',52:'30 Jan',
# #      62:'4 Feb',78:'9 Feb',82:'14 Feb',92:'18 Feb'}
# mm = {2:'5 Jan',14:'11 Jan',42:'25 Jan',62:'4 Feb',84:'15 Feb',104:'25 Feb',124:'6 Mar'}
# for d in mm:
#     path = TextPath((5,0),mm[d])
#     plt.plot(trac['lons'][d],trac['lats'][d],marker='D',markersize=6,color='k')
#     plt.plot(trac['lons'][d],trac['lats'][d],marker=path,markersize=100,color='red')
# if figsav:
#     plt.savefig(join('figs','trajectory_18Mar.png'),**figargs)
# plt.show()