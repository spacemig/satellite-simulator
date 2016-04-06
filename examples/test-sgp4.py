# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 11:22:44 2013

@author: miguelnunes
"""

from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv

line1 = ('1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753')
line2 = ('2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667')

#ISS
#line1 = ('1 25544U 98067A   08264.51782528 −.00002182  00000-0 -11606-4 0  2927')
#line2 = ('2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537')


#HiakaSat tentative
#line1 = ('1 25544U 98067A   08264.51782528 −.00002182  00000-0 -11606-4 0  2927')
#line2 = ('2 25544  91.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537')

satellite = twoline2rv(line1, line2, wgs72)

position, velocity = satellite.propagate(2000, 6, 29, 12, 50, 19)