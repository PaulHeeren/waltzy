# -*- coding: utf-8 -*-
"""
CLEAN THE WAVCAL-FILES FOR THE CERES PIPELINE.
The files were created before with create_wavdats.py.
They still include columns for the amplitudes of the lines,
their FWHM, and their approximated wavelengths.
In the end, they should only hold the number of lines (first col),
the pixel positions within the order (as int, second col),
the exact wavelength (in Angstrom, third col),
and the type of line (fourth col).

Created on Mon Oct  9 10:22:31 2017

@author: pheeren
"""

import glob
import sys
import numpy as np

input_dir   = "wavcals_init/"   # here are the input wavdats
output_dir  = "wavcals/"        # here are the final output wavdats

# read all wavdats into a list
wavdats = glob.glob(input_dir + "*iwdat")

print "Reading data from files and creating new wavcals:"

for dat in wavdats:
    print dat
    readfile = open(dat, 'r')
    # read the data into an array
    alllines = []
    for line in readfile:    
        linelist = [elt.strip() for elt in line.split()]
        # most lines have 7 columns (no. of lines, pixel position, amplitude,
        # width, exact wavelength, approx. wavelength and type of line),
        # but some have 6 (no exact wavelength).
        # However, only the ones with exact wavelength are kept.
        if len(linelist) == 7:
            alllines.append(linelist)    
    
    readfile.close()
    
    alllines = np.array(alllines)
    #print alllines.shape
    #print alllines
    
    #alllines[:,1] = float(alllines[:,1])
    #alllines[:,0] = int(alllines[:,0])
    #alllines[:,1] = np.rint(alllines[:,1])
    
    writefile = open(output_dir + dat[12:], 'w')
    for line in range(alllines.shape[0]):
        writefile.write(alllines[line,0] + "\t" + str(int(np.rint(float(alllines[line,1])))) + "\t" + \
                    alllines[line,4] + "\t" + alllines[line,6] + "\n")
    writefile.close()