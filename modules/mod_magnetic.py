# make sure you have geomag installed, if not then there are two ways to deal with it
# 1) $python setup.py install for first time
# 2) import the module from the thirdparty 
import geomag 

from math import degrees 
from numpy import sqrt
from mod_frames import ecef_from_ned

# magnetic field in 
# before was mag_EFEC
def b_ecef(lat,lon,alt):
    # lat, lon in radians
    geo = geomag.geomag.GeoMag()
    #lat = lat_deg*pi/180. #rad 
    #lon = lon_deg*pi/180. #rad
    #alt = alt_km  #km, use km as default

    #deg, km
    # convert rad to deg
    # check using http://www.ngdc.noaa.gov/geomag-web/#igrfwmm
    mag = geo.GeoMag(degrees(lat),degrees(lon),alt) #*180./pi
    
    bx = mag.bx #N
    by = mag.by #E
    bz = mag.bz #D
    
    #print "lat:{} deg".format(mag.lat)
    #print "lon:{} deg".format(mag.lon)
    #print "alt:{} km".format(mag.alt)
    #print "bh:{} nT".format(mag.bh)
    #print "N (bx):{} nT".format(bx)
    #print "E (by):{} nT".format(by)
    #print "D (bz):{} nT".format(bz)
    #print "bz:{} nT".format(mag.bz)
    
    total_b = sqrt(bx**2 + by**2 + bz**2)
    #print "total:{} nT".format(total_b)
    
    
    # python
    #xx = cos(lat*pi/180.) * cos(lon*pi/180.) * (r_earth + alt)
    #yy = cos(lat*pi/180.) * sin(lon*pi/180.) * (r_earth + alt)
    #zz = sin(lat*pi/180.) * (r_earth + alt) # z is 'up'
        #print xx,yy,zz
    #xx = cos(phi)*cos(theta)*r
    #print xx
    
    
    # convert form NED to ECEF
    # normalize
    n = bx/total_b
    e = by/total_b
    d = bz/total_b
    bx_ecef, by_ecef, bz_ecef = ecef_from_ned(lat,lon,alt,n,e,d)
    
    #return the total magnetic field in Tesla
    return bx_ecef*float(total_b)*1e-9, by_ecef*float(total_b)*1e-9, bz_ecef*float(total_b)*1e-9, total_b*1e-9
    