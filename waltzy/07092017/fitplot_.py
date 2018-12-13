# -*- coding: utf-8 -*-
"""
Plot column or row of .fit file or just get information

Created on Wed Aug 30 12:09:19 2017

@author: pheeren
"""

import pyfits
import argparse

import matplotlib
import numpy as np
#matplotlib.use("Agg") 
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#from matplotlib.ticker import NullFormatter  # useful for `logit` scale


# receive input parameters from terminal

parser = argparse.ArgumentParser()
parser.add_argument('fitsfile')
parser.add_argument('-row', default='-1')
parser.add_argument('-col', default='-1')
parser.add_argument('-log', action="store_true", default=False)
parser.add_argument('-s', action="store_true", default=False)
parser.add_argument('-header', action="store_true", default=False)
#parser.add_argument('-axes', default='01')

args = parser.parse_args()
fits = args.fitsfile            # input file
row  = int(args.row)            # row to plot
col  = int(args.col)            # column to plot
log  = args.log                 # logarithmic scale?
save = args.s                   # save picture as .svg?
showheader = args.header        # show header of file
#whichaxes  = args.axes          # which axes contain the information

h = pyfits.open(fits)
firstdata = h[0].data

# This is only for fits-files with 3 dimensions
fitsdata = firstdata[0,:,:]

if showheader == True:
    fitsheader = h[0].header
    print fitsheader

# First case: No col or row specified, give image info and plot matrix
if row == -1 and col == -1:
    h.info()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if log == True:
        #fitsdata[np.where(fitsdata < 0)] = 0.1
        cax = ax.matshow(fitsdata, norm=LogNorm(vmin=fitsdata.min(), vmax=fitsdata.max()))
    else: cax = ax.matshow(fitsdata)
    
    plt.ylabel('Rows')
    plt.xlabel('Columns')
    ax.xaxis.set_ticks_position('bottom')
    ax.invert_yaxis()
    plt.title(fits)
    fig.colorbar(cax)
    
    if save == True:
        name = fits + ".svg"
        plt.savefig(name, dpi=1200)
    
    plt.show()


# Second case: Plot col
elif row == -1 and col != -1:
    plotdata = fitsdata[:,col,0]
    x = fitsdata.shape[0]
    plt.plot(range(x), plotdata)
    if log == True:
        plt.yscale('log')
        
    plt.ylabel('Counts')
    plt.xlabel('Row pixels')
    title = 'Column No. ' + str(col) + ", " + fits
    plt.title(title)
    
    if save == True:
        name = fits + "_col" + str(col) + ".svg"
        plt.savefig(name, dpi=1200)
        
    plt.show()


# Third case: Plot row
elif col == -1 and row != -1:
    plotdata = fitsdata[row,:]
    x = fitsdata.shape[1]
    plt.plot(range(x), plotdata)
    if log == True:
        plt.yscale('log')
        
    plt.ylabel('Counts')
    plt.xlabel('Column pixels')
    title = 'Row No. ' + str(row) + ", " + fits
    plt.title(title)
    
    if save == True:
        name = fits + "_row" + str(row) + ".svg"
        plt.savefig(name, dpi=1200)
        
    plt.show()


# Fourth case: Give counts for one pixel
elif col != -1 and row != -1:
    pix = "(" + str(row) + "," + str(col) + ")"
    print "Count in pixel " + pix + ": ", fitsdata[row, col]