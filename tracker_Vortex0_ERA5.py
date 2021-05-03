#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajectory of Canadian vortex 0
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
pts =   list(np.arange(330,500,5))
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
with gzip.open('ERA5-extract-0-PT5.pkl','rb') as f:
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
lat = 62
lon = -80
upper = 30
lower = 12
lower = 2
upper= 4
lat = 72
lon = -115
date = datetime(2017,8,14,12)
ioff = 2
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
        upper = min(33,trac['kz'][-1] + 3)
        lower = trac['kz'][-1] - 2
    except: pass

    if (var == 'LPV'):
        if i >= 17:
           lower = max(lower,11)
        elif i>= 12:
           lower = max(lower,9)
        elif i>= 10:
           lower = max(lower,7)
        elif i>= 8:
           lower = max(lower,5)
        elif i>= 6:
           lower = max(lower,3)
    if (var == 'O3prime'):
        if i == 11:
            lat = 62
            lon = -91
        if i == 53:
            lat = 48
            lon = 5
            idx = 2
        if i >= 17:
           lower = max(lower,11)
        elif i>= 12:
           lower = max(lower,9)
        elif i>= 10:
           lower = max(lower,7)
        elif i>= 8:
           lower = max(lower,5)
        elif i>= 6:
           lower = max(lower,3)

    #     if i==1:
    #         lon = -78
    #         lat = 60

    dats[i1].var['LPV'] = dats[i1].var['PV']*LaitFactor[:,None,None]
    meanO3 = np.mean(dats[i1].var['O3'],axis=(1,2))
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
        trac_PV[k] = trac[k][:174].copy()
    trac1 = trac_PV
elif var == 'O3prime':
    trac_O3 = {}
    for k in trac.keys():
        trac_O3[k] = trac[k][:174].copy()
    trac2 = trac_O3

#%% print the positions as a function of time
for i in range(len(trac1['dates'])):
    # kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    print(i,trac1['dates'][i],trac1['lons'][i],trac1['lats'][i],'{:2.1f}'.format(trac1['z'][i]),trac1['kz'][i])
    print(i,trac2['dates'][i],trac2['lons'][i],trac2['lats'][i],'{:2.1f}'.format(trac2['z'][i]),trac2['kz'][i])

pickle.dump([trac_PV,trac_O3],open('Vortex0-track.pkl','wb'))

#%% Total displacement
[trac1,trac2] = pickle.load(open('Vortex0-track.pkl','rb'))
# with PV
lats = np.array(trac1['lats'])
lons = np.array(trac1['lons'])
dy = lats[1:]-lats[:-1]
dx = (lons[1:]-lons[:-1])*np.cos(np.deg2rad(0.5*(lats[1:]+lats[:-1])))
ds = np.sqrt(dx**2+dy**2)
print('total path PV ',(2*np.pi*6371/360)*np.sum(ds))
# with O3
lats = np.array(trac2['lats'])
lons = np.array(trac2['lons'])
dy = lats[1:]-lats[:-1]
dx = (lons[1:]-lons[:-1])*np.cos(np.deg2rad(0.5*(lats[1:]+lats[:-1])))
ds = np.sqrt(dx**2+dy**2)
print('total path O3 ',(2*np.pi*6371/360)*np.sum(ds))

#%% Checker

#for i in range(50,len(trac1['dates'])):
images = []
movieOn = True
for i in range(0,82):
    # generates scaled Montogomery potential
    Mo = (cst.g * dats[i+ioff].var['Z'] + cst.Cp * dats[i+ioff].var['T']) * 1.e-3
    kz1 = trac1['kz'][i]
    kz2 = trac2['kz'][i]
    fig = plt.figure(figsize=(11,8))
    projplate = ccrs.PlateCarree(central_longitude=0)
    gs = fig.add_gridspec(2,1)
    ax1 = fig.add_subplot(gs[0,0],projection=projplate)
    ax1 = dats[i+ioff].show('PV',kz1,figsize=None,axf=ax1,show=False,scale=1.e6,xylim=True,
                 txt='PV (PVU)  '+str(i)+'  '+str(kz1)+trac1['dates'][i].strftime('  %Y %b %d %HUTC  ')+str(pts[kz1])+'K'+
                 '   lon'+str(trac1['lons'][i])+'   lat'+str(trac1['lats'][i]))
    ax1.plot(trac1['lons'][i],trac1['lats'][i],'black',marker='x',ms=21,mew=2)
    Mof = gaussian_filter(Mo[kz1,...],4)
    Mof = Mof - np.mean(Mof)
    cs = ax1.contour(Mof,transform=projplate,extent=(dats[i+ioff].attr['lons'][0],
                dats[i+ioff].attr['lons'][-1],dats[i+ioff].attr['lats'][0],dats[i+ioff].attr['lats'][-1]),origin='lower')
    plt.clabel(cs)
    ax2 = fig.add_subplot(gs[1,0],projection=projplate)
    ax2 = dats[i+ioff].show('O3',kz2,figsize=None,axf=ax2,show=False,scale=1.e6,xylim=True,
                 txt='O3 (ppmv)  '+str(i)+'  '+str(kz2)+trac2['dates'][i].strftime('  %Y %b %d %HUTC  ')+str(pts[kz2])+'K'+
                 '   lon'+str(trac2['lons'][i])+'   lat'+str(trac2['lats'][i]))
    Mof = gaussian_filter(Mo[kz2,...],4)
    Mof = Mof - np.mean(Mof)
    cs = ax2.contour(Mof,transform=projplate,extent=(dats[i+ioff].attr['lons'][0],
                dats[i+ioff].attr['lons'][-1],dats[i+ioff].attr['lats'][0],dats[i+ioff].attr['lats'][-1]),origin='lower')
    plt.clabel(cs)
    ax2.plot(trac2['lons'][i],trac2['lats'][i],'black',marker='x',ms=21,mew=2)
    ax1.text(-0.05,1.05,'a)',transform=ax1.transAxes,fontsize=22)
    ax2.text(-0.05,1.05,'b)',transform=ax2.transAxes,fontsize=22)
    if movieOn:
        plt.savefig('tmpim.png',bbox_inches='tight')
        images.append(imageio.imread('tmpim.png'))
    plt.show()
#%%
if movieOn:
    imageio.mimsave('movie_Vortex0.gif', images, fps=2)
    # reshape images to have the same size (required by mp4), the problem occurs in lon
    # find the min size
    # this can certainly be improved to be more pythonesq
    minX = 20000
    for j in range(len(images)):
        minX = min(minX,images[j].shape[1])
    for j in range(len(images)):
        if images[j].shape[1] > minX:
            images[j] = images[j][:,:minX,:]

    imageio.mimsave('movie_Vortex0.mp4', images, fps=2)

#%% Maxichecker
images = []
movieOn = True
# Generate uniform PT field and lats & lons grids
PT = dats[0].var['T'].copy()
for l in range(len(dats[0].attr['levs'])):
    PT[l,...] = dats[0].attr['levs'][l]
for i in range(0,82):
#for i in [12,24]:
    # generates scaled Montogomery potential
    Mo = (cst.g * dats[i+ioff].var['Z'] + cst.Cp * dats[i+ioff].var['T']) * 1.e-3
    kz1 = trac1['kz'][i]
    kz2 = trac2['kz'][i]
    lat1 = trac1['lats'][i]
    lat2 = trac2['lats'][i]
    lon1 = trac1['lons'][i]
    lon2 = trac2['lons'][i]
    z1 = trac1['z'][i]
    z2 = trac2['z'][i]
    fig = plt.figure(figsize=(30,11))
    lats = np.tile(dats[i+ioff].attr['lats'],(33-2,1))
    lons = np.tile(dats[i+ioff].attr['lons'],(33-2,1))
    projplate = ccrs.PlateCarree(central_longitude=0)
    gs = fig.add_gridspec(14,10)
    ax1 = fig.add_subplot(gs[:4,:5],projection=projplate)
    ax1 = dats[i+ioff].show('PV',kz1,figsize=None,axf=ax1,show=False,scale=1.e6,xylim=True,
                 txt='PV (PVU)  '+str(i)+'  '+str(kz1)+trac1['dates'][i].strftime('  %Y %b %d %HUTC  ')+str(pts[kz1])+'K'+
                 '   lon'+str(lon1)+'   lat'+str(lat1))
    Mof = gaussian_filter(Mo[kz1,...],4)
    Mof = Mof - np.mean(Mof)
    cs = ax1.contour(Mof,transform=projplate,extent=(dats[i+ioff].attr['lons'][0],
                dats[i+ioff].attr['lons'][-1],dats[i+ioff].attr['lats'][0],dats[i+ioff].attr['lats'][-1]),origin='lower')
    plt.clabel(cs)
    ax1.plot(lon1,lat1,'white',marker='x',ms=21,mew=4)
    ax2 = fig.add_subplot(gs[:4,5:],projection=projplate)
    ax2 = dats[i+ioff].show('O3',kz2,figsize=None,axf=ax2,show=False,scale=1.e6,xylim=True,
                 txt='O3 (ppmv)  '+str(i)+'  '+str(kz2)+trac2['dates'][i].strftime('  %Y %b %d %HUTC  ')+str(pts[kz2])+'K'+
                 '   lon'+str(lon2)+'   lat'+str(lat2))
    Mof = gaussian_filter(Mo[kz2,...],4)
    Mof = Mof - np.mean(Mof)
    cs = ax2.contour(Mof,transform=projplate,extent=(dats[i+ioff].attr['lons'][0],
                dats[i+ioff].attr['lons'][-1],dats[i+ioff].attr['lats'][0],dats[i+ioff].attr['lats'][-1]),origin='lower')
    plt.clabel(cs)
    ax2.plot(lon2,lat2,'white',marker='x',ms=21,mew=4)
    ax3 = fig.add_subplot(gs[5:9,1:5])
    ax3.set_anchor('C',share=True)
    ax3 = dats[i+ioff].chartlatz('PV',lon1,levs=(2,32),show=False,figsize=None,axf=ax3,scale=1.e6,
                 txt='PV (PVU) '+str(i)+' lon '+str(lon1)+trac1['dates'][i].strftime('                                    %Y %b %d %HUTC  ') )
    ax3.plot(lat1,z1,'white',marker='x',ms=21,mew=4)
    idx3 = np.where(dats[i+ioff].attr['lons'] == lon1)[0][0]
    cf3 = ax3.contour(lats,dats[i+ioff].var['Z'][2:33,:,idx3]/1000,PT[2:33,:,idx3],
           colors='white')
    plt.clabel(cf3,cf3.levels[1:-1:2],inline=True, fontsize=16,fmt='%0.f')
    ax4 = fig.add_subplot(gs[5:9,6:])
    ax4 = dats[i+ioff].chartlatz('O3',lon2,levs=(2,32),show=False,figsize=None,axf=ax4,scale=1.e6,
                 txt='O3 (ppmv)  '+str(i)+' lon '+str(lon2)+trac1['dates'][i].strftime('                                 %Y %b %d %HUTC  ') )
    ax4.plot(lat2,z2,'white',marker='x',ms=21,mew=4)
    idx4 = np.where(dats[i+ioff].attr['lons'] == lon2)[0][0]
    cf4 = ax4.contour(lats,dats[i+ioff].var['Z'][2:33,:,idx4]/1000,PT[2:33,:,idx4],
           colors='white')
    plt.clabel(cf4,cf4.levels[1:-1:2],inline=True, fontsize=16,fmt='%0.f')
    ax5 = fig.add_subplot(gs[10:,1:5])
    ax5 = dats[i+ioff].chartlonz('PV',lat1,levs=(2,32),show=False,figsize=None,axf=ax5,scale=1.e6,
                 txt='PV (PVU) '+str(i)+' lat '+str(lat1)+trac1['dates'][i].strftime('                                       %Y %b %d %HUTC  ') )
    ax5.plot(lon1,z1,'white',marker='x',ms=21,mew=4)
    idx5 = np.where(dats[i+ioff].attr['lats'] == lat1)[0][0]
    cf5 = ax5.contour(lons,dats[i+ioff].var['Z'][2:33,idx5,:]/1000,PT[2:33,idx5,:],
           colors='white')
    plt.clabel(cf5,cf5.levels[1:-1:2],inline=True, fontsize=16,fmt='%0.f')
    ax6 = fig.add_subplot(gs[10:,6:])
    ax6 = dats[i+ioff].chartlonz('O3',lat2,levs=(2,32),show=False,figsize=None,axf=ax6,scale=1.e6,
                 txt='O3 (ppmv) '+str(i)+' lat '+str(lat2)+trac1['dates'][i].strftime('                                     %Y %b %d %HUTC  ') )
    ax6.plot(lon2,z2,'white',marker='x',ms=21,mew=4)
    idx6 = np.where(dats[i+ioff].attr['lats'] == lat2)[0][0]
    cf6 = ax6.contour(lons,dats[i+ioff].var['Z'][2:33,idx6,:]/1000,PT[2:33,idx6,:],
           colors='white')
    plt.clabel(cf6,cf6.levels[1:-1:2],inline=True, fontsize=16,fmt='%0.f')
    ax1.text(-0.05,1.05,'a)',transform=ax1.transAxes,fontsize=22)
    ax2.text(-0.05,1.05,'b)',transform=ax2.transAxes,fontsize=22)
    ax3.text(-0.15,1.05,'c)',transform=ax3.transAxes,fontsize=22)
    ax4.text(-0.15,1.05,'d)',transform=ax4.transAxes,fontsize=22)
    ax5.text(-0.15,1.05,'e)',transform=ax5.transAxes,fontsize=22)
    ax6.text(-0.15,1.05,'f)',transform=ax6.transAxes,fontsize=22)
    ax3.set_ylim(11,20)
    ax4.set_ylim(11,20)
    ax5.set_ylim(11,20)
    ax6.set_ylim(11,20)
    if movieOn:
        plt.savefig('tmpim.png',bbox_inches='tight')
        images.append(imageio.imread('tmpim.png'))
    plt.show()
#%%
if movieOn:
    imageio.mimsave('movieMax_Vortex0.gif', images, fps=2)
    # reshape images to have the same size (required by mp4), the problem occurs in lon
    # find the min size
    # this can certainly be improved to be more pythonesq
    minX = 20000
    for j in range(len(images)):
        minX = min(minX,images[j].shape[1])
    for j in range(len(images)):
        if images[j].shape[1] > minX:
            images[j] = images[j][:,:minX,:]

    imageio.mimsave('movieMax_Vortex0.mp4', images, fps=2)


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
[trac1,trac2] = pickle.load(open('Vortex0-track.pkl','rb'))
#%%
fig,((ax0,ax1),(ax2,ax3),(ax4,ax5))=plt.subplots(3,2,figsize=(10,10),sharex=True)
pvl1 = np.array(trac1['pvano'])*(np.array(trac1['pt'])/500)**(-4.5)
pvl2 = np.array(trac2['pvano'])*(np.array(trac2['pt'])/500)**(-4.5)
ax0.plot(
         trac2['dates'],trac2['lats'],
         trac1['dates'],trac1['lats'],linewidth=2)
ax1.plot(
         trac2['dates'],np.array(trac2['lons']),
         trac1['dates'],np.array(trac1['lons']),linewidth=2)
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
fig.suptitle('Canadian event   vortex 0 tracking')
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('VortexMotion_Vortex0_ERA5PV.png'),**figargs)
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