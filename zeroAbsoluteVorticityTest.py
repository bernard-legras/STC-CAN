# -*- coding: utf-8 -*-
"""
Test the hypothesis that the vorticity is near -f

Created on Sat Oct 10 10:56:03 2020

@author: Bernard Legras
"""
import pickle
import constants as cst
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10,10))
plt.subplot(2,2,1)
[trac_PV,trac_O3]=pickle.load(open('Vortex0-track.pkl','rb'))
fcorio = 2 * cst.Omega * np.sin(np.deg2rad(trac_PV['lats']))
vratio = 1+trac_PV['vo']/fcorio
plt.plot(trac_PV['dates'],vratio)
plt.ylabel(u'vortex 0 1+$\zeta$/f')
plt.setp(plt.xticks()[1], rotation=30, ha='right')
plt.subplot(2,2,2)
[trac_PV,trac_O3]=pickle.load(open('VortexA-track.pkl','rb'))
fcorio = 2 * cst.Omega * np.sin(np.deg2rad(trac_PV['lats']))
vratio = 1+trac_PV['vo']/fcorio
plt.plot(trac_PV['dates'],vratio)
plt.ylabel(u'vortex A 1+$\zeta$/f')
plt.setp(plt.xticks()[1], rotation=30, ha='right')
plt.subplot(2,2,3)
[trac_PV,trac_O3]=pickle.load(open('VortexB1-track.pkl','rb'))
fcorio = 2 * cst.Omega * np.sin(np.deg2rad(trac_PV['lats']))
vratio = 1+trac_PV['vo']/fcorio
plt.plot(trac_PV['dates'],vratio)
plt.ylabel(u'vortex B1 1+$\zeta$/f')
plt.setp(plt.xticks()[1], rotation=30, ha='right')
plt.subplot(2,2,4)
[trac_PV,trac_O3]=pickle.load(open('VortexB2-track.pkl','rb'))
fcorio = 2 * cst.Omega * np.sin(np.deg2rad(trac_PV['lats']))
vratio = 1+trac_PV['vo']/fcorio
plt.plot(trac_PV['dates'],vratio)
plt.ylabel(u'vortex B2 1+$\zeta$/f')
plt.setp(plt.xticks()[1], rotation=30, ha='right')
#fig.suptitle('Australia 2020')
plt.show()

