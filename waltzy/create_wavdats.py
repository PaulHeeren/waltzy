# -*- coding: utf-8 -*-
"""
CREATE THE WAVCALS-FILES FOR THE CERES PIPELINE.
First take the txt-file with all the lines created by ext_tharanalyze.py
and fit a rough wavelength approximation (based on the spectrograph
equation). The result is shown in comparison to known line wavelengths
put into the file beforehand. If wished, a low order polynomial can be fitted
to achieve a better result. The outcome for each order is then written to 
different order_XXX.iwdat files.

Created on Tue Sep 26 16:11:37 2017

@author: pheeren
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.optimize as optimization

# some parameters
order1 = 54         # first physical order in file
cpix   = 1024.      # central pixel of the orders
gamma  = 0.35       # gamma angle of WALTZY
grat_d = .00000316  # lines per Angstrom of the grating
grat_i = 63.9       # blaze angle of the grating
cam_f  = 530000.    # focal length of camera in mum
ccd_s  = 13.5       # CCD pixel size in mum

degtorad = np.pi/180.   # convert degrees to rad
# directory to write the wavcals to
dirout = "wavcals_init/"

# define a function for (optional) fitting of wavelength approx
def fitfunc(x, a, b):
    return a + b*x


# first take the given filename
if len(sys.argv) == 1:
    input_filename = raw_input('Please enter filename: ')
else:
    input_filename = sys.argv[1]

# all the lines are in a txt-file (and maybe some first wavelengths added by hand)
alltxt = open(input_filename, 'r')

print "Reading data from file."
# read the data into an array
alllines = []
for line in alltxt:    
    linelist = [elt.strip() for elt in line.split()]
    # most lines have only 4 columns (order number, pixel position, amplitude, width),
    # but some have 5 (also a wavelength added by hand).
    # However, we need 6 columns for all of them, so fill the rest up with zeros.
    if len(linelist) == 5:
        linelist.append(0)
    elif len(linelist) == 4:
        linelist.append(0)
        linelist.append(0)
    alllines.append(linelist)
alltxt.close()
alllines = np.array(alllines)

alllines = alllines.astype(np.float)
#print alllines.shape
#print alllines[42]


# now convert the code order (starting at 0) to the physical order (starting at order1)
alllines[:,0] = alllines[:,0] + order1

print "Calculating wavelengths for each line, according to wavelength equation of WALTZY..."
# calculate the approx wavelengths for each order, according to the wavelength equation   
equ_factor = np.cos(gamma*degtorad) / (grat_d)

wav_fixed = []
wav_approx = []

for line in range(alllines.shape[0]):
    # offset of line from central pixel position
    dpix = cpix - alllines[line,1]
    # outgoing angle from echelle grating for the line (f of camera: 530mm)
    beta = (grat_i*degtorad) + (dpix*ccd_s/cam_f)
    # approx wavelength at pixel position (this is the grating equation)
    alllines[line,5] = equ_factor * (np.sin(grat_i*degtorad) + np.sin(beta)) / alllines[line,0]
    
    # lines that already have an associated wavelength from before are used for comparison
    if alllines[line,4] != 0:
        wav_fixed.append(alllines[line,4])
        wav_approx.append(alllines[line,5])

wav_fixed = np.array(wav_fixed)
wav_approx = np.array(wav_approx)
plt.scatter(wav_fixed, (wav_approx-wav_fixed))
plt.axhline(y=0, linewidth=2, color='r')
plt.show()

# User wishes fitting procedure for better results?
fitting_wish = raw_input('You want a simple fit? (Y for Yes, N for No) ')

# if fitting wished, do it and show the results
par0 = np.array([0., 1.])
if fitting_wish == 'Y':
    print "Yes Master, fitting wavelengths..."
    par, errors = optimization.curve_fit(fitfunc, wav_approx, wav_fixed, par0)
    print "Fit output parameters (y = a + x*b): ", par
    
    print "Using these parameters on the wavelength approximation..."
    alllines[:,5] = par[0] + alllines[:,5] * par[1]
    wav_approx = par[0] + wav_approx * par[1]
    plt.scatter(wav_fixed, (wav_approx-wav_fixed))
    plt.axhline(y=0, linewidth=2, color='r')
    plt.show()
    
elif fitting_wish == "N":
    print "Aye aye, no fitting..."

print "Write the wavcal files into the directory: ", dirout
# now write to wavcal files
name = []
order = order1 - 1  # running variable for file creation

for line in range(alllines.shape[0]):
    if int(alllines[line,0]) != order:
        order += 1
        if order < 100:
            name = "order_0" + str(order) + ".iwdat"
        else:
            name = "order_" + str(order) + ".iwdat"
        outfile = open(dirout+name, 'w')
    else:
        outfile = open(dirout+name, 'a')
    
    outfile.write("1")
    for a in range(alllines.shape[1]-1):
        outfile.write("\t" + str(alllines[line,a+1]))
    outfile.write("\n")
    outfile.close()