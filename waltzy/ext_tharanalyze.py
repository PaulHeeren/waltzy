# -*- coding: utf-8 -*-
"""
Analyse Thorium-Argon-Specs for widths and positions of their lines.
Need traced and extracted spectrum for this.
The lines are then written to a txt file, holding the following parameters:
- order number of the line (as in code, not physical)
- fitted pixel position of the line within the order
- fitted amplitude of the line
- fitted width of the line (in pixels, as FWHM)

Created on Fri Sep 15 10:06:02 2017

@author: pheeren
"""

import astropy.io.fits as pyfits
import numpy as np
import scipy
from scipy import optimize

import matplotlib.pyplot as plt

import pickle
import sys

# define some parameters
linestrength = 200  # min amplitude for lines; everything below is not accepted


# define the fitfunction (gaussian)
def fitfunc(p,x):
    ret = p[0] + p[1] * np.exp(-.5*((x-p[2])/p[3])**2)
    return ret
errfunc = lambda p,y,x: np.ravel( (fitfunc(p,x)-y) )


# start with spectrum analysis: first take the given filename
if len(sys.argv) == 1:
    input_filename = raw_input('Please enter filename:')
else:
    input_filename = sys.argv[1]

# Open spec to analyse
h = pyfits.open(input_filename)
tharspec = h[0].data[:,1,:]     # this is the format of the fits-files produced in CERES
print "File found, starting..."

print tharspec.shape
sys.exit

amps = []       # fitted amplitudes of lines
refs = []       # fitted positions of lines
wids = []       # fitted widths of lines
order = []      # order numbers of lines
nr_of_lines = 0
# Loop through all the orders
for i in range(tharspec.shape[0]):
    
    thardata = tharspec[i,:]
    print "Analyse order: ", i
    
    refw = 8        # region around each pixel to analyse to find emission lines
    sigc = 1.5      # 1st width of gaussian to find objects in the spectra
        
    ejx = np.arange(len(thardata))
    ccf=[]
    j = 0
    # Now loop over all pixels to find possible emission lines
    while j < len(thardata):
        # pixels too close to the edges of the orders are not accepted
        if j - refw < 0:
            ccf.append(0.)
        elif j + refw + 1 > len(thardata):
            ccf.append(0.)
        else:
            x = ejx[j-refw:j+refw+1]
            y = thardata[j-refw:j+refw+1]
            # fit gaussian
            g = np.exp(-0.5*((x-j)/sigc)**2)
            ccf.append(np.add.reduce(y*g))
        j+=1
    
    j = 1
    maxs = []
    while j < len(ccf)-2:
        # only take pixels that have highest flux - these are peaks!
        if ccf[j]>ccf[j-1] and ccf[j]>ccf[j+1]:
            maxs.append(j)
        j+=1
    
    maxs = np.array(maxs)
    ccf = np.array(ccf)
    pos = np.arange(len(ccf))[maxs]

    # only lines with high emission above a certain threshold are accepted
    I = np.where(thardata[pos] > linestrength)[0]
    
    pos = pos[I]
    nr_of_lines += len(pos)
        
    #print pos    
    """plt.plot(thardata)
    plt.plot(ccf)
    plt.plot(pos,thardata[pos],'r.')
    plt.show()"""
    
    # now loop throuh all the lines in the order and fit gaussians
    exap2 = 4         # region around each identified line to fit
    dev = 1.5         # input width of gaussian
    for j in range(len(pos)):
        x = ejx[int(pos[j]-exap2):int(pos[j]+exap2+1)]
        y = thardata[int(pos[j]-exap2):int(pos[j]+exap2+1)]
        
        tx1 = np.arange(x[0]-dev,x[0],1)
        tx2 = np.arange(x[-1]+1,x[-1]+dev+1,1)
        ty1 = np.zeros(len(tx1))
        ty2 = np.zeros(len(tx2))
        x = np.hstack((tx1,x,tx2))
        y = np.hstack((ty1,y,ty2))
        y -= y.min()
        
        # fit gaussian
        p, success =  scipy.optimize.leastsq(errfunc, [y.min(),y.max()-y.min(),x.mean(),dev], args=(y,x))
        # append fitting parameters
        amps.append(p[1])
        refs.append(p[2])
        wids.append(p[3])
        order.append(i)

amps = np.array(amps)
refs = np.array(refs)
wids = np.array(wids)
order = np.array(order)
wids = 2.355*wids     # convert to fwhm

# calculate resolution
factor = 2*np.tan(63.9*np.pi/180)*530000
resolution = np.array(factor/(wids*13.5))

print 'Total number of identified lines: ', nr_of_lines
avgres = np.sum(resolution)/nr_of_lines
print 'Average resolution: ', avgres

# create one array that holds all the parameters
alllines = []
for line in range(nr_of_lines):
    # only lines that are not overflowing
    if wids[line] < 6.:
        alllines.append([order[line], refs[line], wids[line], resolution[line], amps[line]])

alllines = np.array(alllines)
# save array data for possible later use
pickle.dump(alllines, open(input_filename.split('/')[-1][:-7]+'_lines.pkl', 'w'))

# Write the lines to a txt-file (not the resolution)
txtfile = open(input_filename.split('/')[-1][:-7]+"_lines.txt", 'w')
for line in range(alllines.shape[0]):
    txtfile.write("{}\t{}\t{}\t{:.2f}\n".format(int(alllines[line,0]), alllines[line,1], \
                                            int(alllines[line,4]), alllines[line,2]))
txtfile.close()

# plot results
plt.figure()
plt.scatter(alllines[:,1], alllines[:,0], s = 50, c = alllines[:,2], alpha = 0.7, vmin = 0, vmax = 6, cmap = 'hot')
plt.colorbar()
plt.title('Testing focus - FWHM')
plt.show()

plt.hist(alllines[:,2], bins=18)
plt.xlabel('Line width')
plt.title('Testing focus - FWHM')
plt.show()