# -*- coding: utf-8 -*-
# Python Space Simulator
# Test IGRF

# import system modules
import os
import sys
from math import degrees, radians, pi, sin, cos, sqrt

# get the current script file path
this_folder     = os.path.dirname(os.path.abspath(__file__))+'/'
root_folder_pss = os.path.dirname(os.path.abspath(__file__+'./..'))+'/'
modules_folder  = root_folder_pss+'modules/'
geomag_folder   = root_folder_pss+'thirdparty/geomag-0.9/'

# add the folders to the system path 
if modules_folder not in sys.path:
    sys.path.append(modules_folder)
    
if geomag_folder not in sys.path:
    sys.path.append(geomag_folder)
    
from mod_magnetic import b_ecef #orbital dynamics support library
from mod_frames import ecef_from_geodetic
import geomag
from datetime import date

lat, lon, alt = 0, 0, 0
time = date(2000,1,1)

geo = geomag.geomag.GeoMag()
mag = geo.GeoMag(degrees(lat), 
                degrees(lon),
                alt,
                time) #*180./pi
bx, by, bz = mag.bx, mag.by, mag.bz #N, #E, #D

bx_ecef, by_ecef, bz_ecef, b_total = b_ecef(lat,lon,alt)
xx_ecef, yy_ecef, zz_ecef = ecef_from_geodetic(lat,lon,alt)

print bx, by, bz
print bx_ecef, by_ecef, bz_ecef, b_total
print xx_ecef, yy_ecef, zz_ecef


# compare to matlab code by Charles L. Rino 
# xyzt=igrf11syn(2000,0,0,0)

bx_matlab = 2.746494738862072e4
by_matlab = -0.350415289549645e4
bz_matlab = -1.482775854948421e4

bx_matlab = 2.746494559029745e4
by_matlab = -0.350415320781053e4 
bz_matlab = -1.482776137319431e4

print "Error with Matlab: {:.1f} %, {:.1f} %, {:.1f} %".format(
abs(1-bx_matlab/bx)*100,
abs(1-by_matlab/by)*100,
abs(1-bz_matlab/bz)*100)
