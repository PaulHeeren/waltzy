# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 16:31:41 2017

@author: pheeren
"""

import pyfits
import argparse

import matplotlib
import numpy as np
#matplotlib.use("Agg") 
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ImgList = np.array(["continuum-15s.fit", "led-0010.fit"])

n = len(ImgList)
h = pyfits.open(ImgList[0])
d1 = h[0].data

if n > 1:
    for i in range(n-1):
        td = pyfits.open(ImgList[i+1])
        d1 = np.dstack((d1,td[0].data))

    d1 = np.median(d1,axis=2)
    
fitsdata = d1

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(fitsdata)
    
plt.ylabel('Rows')
plt.xlabel('Columns')
ax.xaxis.set_ticks_position('bottom')
fig.colorbar(cax)
   
plt.show()