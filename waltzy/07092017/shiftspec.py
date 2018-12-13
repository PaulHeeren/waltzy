# -*- coding: utf-8 -*-
"""
Shifts spec by some pixels or so - either in dispersion direction
or in cross-dispersion direction (to test
the line detection algorithm of the CERES pipeline).

Created on Thu Sep 14 17:14:31 2017

@author: pheeren
"""

import pyfits
import numpy as np

h = pyfits.open("continuum_001.fits")
d = h[0].data
d1 = d[0,:,:]

dnew = np.zeros(d.shape)

# Shift spectrum in cross-dispersion direction
for i in range(d1.shape[1]):
    if i < d1.shape[1] - 4:
        dnew[0,:,i] = d1[:,i+4]
    else:
        dnew[0,:,i] = d1[:,i]
"""
# Shift spectrum in dispersion direction
for i in range(d1.shape[0]):
    if i < d1.shape[0] - 50:
        dnew[0,i,:] = d1[i+50,:]
    else:
        dnew[0,i,:] = d1[i,:]"""

hdu = pyfits.PrimaryHDU( dnew )
hdu.writeto( "continuum_001_shiftc.fits" )
    