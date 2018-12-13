# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 16:43:11 2017

@author: pheeren
"""


import sys
from pylab import *

base = '../'
sys.path.append(base+"utils/Continuum")
sys.path.append(base+"utils/Correlation")
sys.path.append(base+"utils/OptExtract")
sys.path.append(base+"utils/GLOBALutils")

baryc_dir= base+'utils/SSEphem/'
sys.path.insert(1, base+"utils/SSEphem/")
ephemeris='DEc403'

import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt

# ecpipe modules
import continuum
import correlation
import Marsh
import waltzyutils
import GLOBALutils

# other useful modules
import pyfits
import pickle
import os
import numpy as np
import scipy
import scipy.interpolate
from math import radians as rad
import argparse
import warnings
warnings.filterwarnings("ignore")

import ephem
import jplephem

from matplotlib.backends.backend_pdf import PdfPages

import statsmodels.api as sm
lowess = sm.nonparametric.lowess

# receive input parameters from terminal

parser = argparse.ArgumentParser()
parser.add_argument('directorio')
parser.add_argument('-o2do',default='all')
parser.add_argument('-just_extract', action="store_true", default=False)
parser.add_argument('-do_class', action="store_true", default=False)
parser.add_argument('-avoid_plot', action="store_true", default=False)
parser.add_argument('-npools', default=1)
parser.add_argument('-reffile',default='default')
parser.add_argument('-dirout',default='default')

args = parser.parse_args()
DoClass     = args.do_class
avoid_plot  = args.avoid_plot
dirin       = args.directorio
object2do   = args.o2do
JustExtract = args.just_extract
npools      = int(args.npools)
reffile     = args.reffile
dirout      = args.dirout

# prepare directories

if dirin[-1] != '/':
    dirin = dirin + '/'

if dirout == 'default':
    dirout = dirin[:-1]+'_red/'

if not os.access(dirout,os.F_OK):
    os.system('mkdir '+dirout)
if os.access(dirout+'proc',os.F_OK):
    os.system('rm -r '+dirout+'proc')
os.system('mkdir '+dirout+'proc')

f_res = open(dirout+'proc/'+'results.txt','w')

if reffile == 'default':
    reffile = dirin+'reffile.txt'
    
    ####### GLOBAL VARIABLES #####
force_pre_process  = False
force_flat_extract = False
force_sci_extract  = False
force_thar_extract = False	
force_tharxc       = False
force_thar_wavcal  = False
force_spectral_file_build = True
force_stellar_pars = False
dumpargon          = False
minlines_glob      = 700

Inverse_m          = True
use_cheby          = True
MRMS               = 100    # max rms in m/s, global wav solution

trace_degree       = 5      # no. of coefficients of polynomial to fit the order traces
Marsh_alg          = 0      # use Marsh (0) or Horne (1) optimal extraction algorithm
ext_aperture       = 8      # width of the echelle orders in pixels
NSigma_Marsh       = 5      # no. of standard dev. to ignore cosmics in weight calculation
NCosmic_Marsh      = 5      # same as above, but for optimal extraction algorithm
S_Marsh            = 0.4    # fraction of a pixel considered for the interpolation
N_Marsh            = 3      # grado polinomio 
min_extract_col    = 0      # min pixel to start order extraction (in dispersion direction)
max_extract_col    = 2044   # max pixel for order extraction

ncoef_x            = 3
ncoef_m            = 8
npar_wsol = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2

# models_path = base+"data/COELHO_MODELS/R_40000b/"  # for atmospheric modelling
order_dir   = base+"waltzy/wavcals/"   # Th-Ar template specs for each order

n_useful = 47               # no. useful orders
ro0      = 68               # lowest physical order for extraction

if True:
    ra          = 213.91529 #360.*float(h[0].header['RA'])/(3600.*24)
    dec         = 19.18242 #float(h[0].header['DEC'])/3600
    altitude    = 560. #float(h[0].header['HIERARCH CAHA TEL GEOELEV'])
    latitude    = 49.399 #float(h[0].header['HIERARCH CAHA TEL GEOLAT'])
    longitude   = 8.721 #float(h[0].header['HIERARCH CAHA TEL GEOLON'])
    epoch       = 2451545. #float(h[0].header['EQUINOX'])
    
    mjd,mjd0 = 57914.895856, 17914.895856

    """ra2,dec2 = GLOBALutils.getcoords(obname, mjd, filen=reffile)
    if ra2 !=0 and dec2 != 0:
        ra = ra2
        dec = dec2
    else:
        print '\t\tUsing the coordinates found in the image header.'
    """

    # Find lambda_bary/lambda_topo using JPLEPHEM
    iers                    = GLOBALutils.JPLiers( baryc_dir, mjd-999.0, mjd+999.0 )
    obsradius, R0           = GLOBALutils.JPLR0( latitude, altitude)
    obpos                   = GLOBALutils.obspos( longitude, obsradius, R0 )
    
    print iers
    print obsradius
    print R0
    print obpos
    
    jplephem.set_ephemeris_dir( baryc_dir , ephemeris )
    jplephem.set_observer_coordinates( float(obpos[0]), float(obpos[1]), float(obpos[2]) )
    
    res         = jplephem.doppler_fraction(float(ra/15.0), float(dec), long(mjd), float(mjd%1), 1, 0.0)
    lbary_ltopo = 1.0 + res['frac'][0]
    bcvel_baryc = ( lbary_ltopo - 1.0 ) * 2.99792458E5
    
    print "\t\tBarycentric velocity:", bcvel_baryc
    
    res  = jplephem.pulse_delay(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
    mbjd = mjd + res['delay'][0] / (3600.0 * 24.0)

    # Moon Phase Calculations
    gobs      = ephem.Observer()  
    gobs.name = h[0].header['TELESCOP']
    gobs.lat  = rad(latitude)  # lat/long in decimal degrees  
    gobs.long = rad(longitude)
    
    gobs.date = h[0].header['DATE'][:10] + ' ' + h[0].header['DATE'][11:]
    mephem    = ephem.Moon()
    mephem.compute(gobs)

    Mcoo = jplephem.object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Mp = jplephem.barycentric_object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Sp = jplephem.barycentric_object_track("Sun", int(mjd), float(mjd%1), 1, 0.0)
    res  = jplephem.object_doppler("Moon", int(mjd), mjd%1, 1, 0.0)
    lunation,moon_state,moonsep,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,ra,dec)
    refvel = bcvel_baryc + moonvel
    print '\t\tRadial Velocity of scattered moonlight: ', refvel
    
    moon_alts.update({fsim:mephem.alt})
    moon_ills.update({fsim:lunation})