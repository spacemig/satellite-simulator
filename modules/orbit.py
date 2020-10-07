#orbital dynamics support library

import numpy as np
from numpy import linalg as la
#from math import pow, degrees 

# $python setup.py install for first time
#import geomag 

# $python setup.py install for first time
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv

from mod_constants import vector_zero

def normalize(a):
    return a/la.norm(a)

def dcos(a,b):
    # compute direction cosine from vectors a and b 
    # then normalize
    return np.dot(a,b)/(la.norm(a)*la.norm(b))
    


def state_satellite(dateTime,tle):
    year, month, day, hour, minute, second = dateTime
    #ISS
    #line1 = ('1 25544U 98067A   08264.51782528  .00002182  00000-0 -11606-4 0  2927')
    #line2 = ('2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537')
    
    #line1 = ('1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753')
    #line2 = ('2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667')
    
    # similar to hiakasat orbit
    #line1 = ('1 39143U 13016B   13112.53874487 -.00003669  11601-4  00000+0 0    49')
    #line2 = ('2 39143 091.6094 306.9444 0018244 348.7202  11.3600 16.10627010   119')
    
    #2 39143  91.6094 306.9444 0018244 348.7202  11.3600 16.10627010   119  
    #2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537
      
    satellite = twoline2rv(tle.line1, tle.line2, wgs72)
    
    #position and velocity are given in ECI True Equator Mean Equinox (TEME) frame
    position, velocity = satellite.propagate(year, month, day, hour, minute, second)
    
    # convert TEME to a standard frame, 
    # suggestion: rotate to Pseudo Earth Fixed (PEF) using 
    # Greenwich Mean Sideral Time (GMST) and then rotate to other standard frame
    # reference: http://www.centerforspace.com/downloads/files/pubs/AIAA-2008-6770.pdf Appendix B
    return position, velocity