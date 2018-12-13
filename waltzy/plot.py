# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 15:04:43 2017

@author: pheeren
"""

import pyfits
import argparse

import matplotlib
import numpy as np
#matplotlib.use("Agg") 
import matplotlib.pyplot as plt
import scipy
import pickle

fits = "continuum-0010.fit"

h = pyfits.open(fits)
fitsdata = h[0].data
fitsdata = np.flip(fitsdata, 0)
fitsdata = np.asarray(zip(*fitsdata[::-1]))
fitsdata = np.asarray(zip(*fitsdata[::-1]))

d = pickle.load(open("trace.pkl", 'r'))
c_all = d['c_all']
nord = d['nord']
p2 = d['p2']     ###################################
Centers = np.zeros((len(c_all), fitsdata.shape[1]))

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(fitsdata)

for i in range(nord):
    Centers[i,:] = scipy.polyval(c_all[i,:],np.arange(len(Centers[i,:])))
    plt.plot(Centers[i,:], np.arange(len(Centers[i,:])), 'k', linewidth=0.4)
    plt.plot(p2[i,:], np.arange(len(Centers[i,:])), 'ro')

plt.ylabel('Rows')
plt.xlabel('Columns')
ax.xaxis.set_ticks_position('bottom')
plt.title(fits)
fig.colorbar(cax)

plt.savefig('continuum-15s-traceplot.svg', dpi=1200)

plt.show()