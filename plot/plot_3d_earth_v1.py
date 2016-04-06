import numpy as np
from math import pow, degrees, radians, pi, sin, cos, sqrt
from scipy import mat,  arctan, arctan2 #pi, cos, sin,

# python setup.py install for first time
import geomag 

ocean_blue = (0.4, 0.5, 1.0)
r_km = 6371.0 # km 
r = r_km

###############################################################################
from mayavi import mlab
mlab.figure(1, bgcolor=(0.48, 0.48, 0.48), fgcolor=(0, 0, 0),
               size=(400, 400))
mlab.clf()

###############################################################################
# Display continents outline, using the VTK Builtin surface 'Earth'
from mayavi.sources.builtin_surface import BuiltinSurface
continents_src = BuiltinSurface(source='earth', name='Continents')
# The on_ratio of the Earth source controls the level of detail of the
# continents outline.
continents_src.data_source.on_ratio = 2
continents_src.data_source.radius = r
continents = mlab.pipeline.surface(continents_src, color=(0, 0, 0))



################################################################################
## Display a semi-transparent sphere, for the surface of the Earth
#
## We use a sphere Glyph, throught the points3d mlab function, rather than
## building the mesh ourselves, because it gives a better transparent
## rendering.
#sphere = mlab.points3d(0, 0, 0, scale_mode='none',
#                                scale_factor=r*2,
#                                color=ocean_blue,
#                                resolution=50,
#                                opacity=0.99,
#                                name='Earth')
#
## These parameters, as well as the color, where tweaked through the GUI,
## with the record mode to produce lines of code usable in a script.
#sphere.actor.property.specular = 0.45
#sphere.actor.property.specular_power = 5
## Backface culling is necessary for more a beautiful transparent
## rendering.
#sphere.actor.property.backface_culling = True

###############################################################################
# Plot the equator and the tropiques
theta = np.linspace(0, 2*np.pi, 100)
for angle in (-np.pi/6, 0, np.pi/6):
    x = r_km*np.cos(theta)*np.cos(angle)
    y = r_km*np.sin(theta)*np.cos(angle)
    z = r_km*np.ones_like(theta)*np.sin(angle)

    mlab.plot3d(x, y, z, color=(1, 1, 1),
                        opacity=0.1, tube_radius=None)



#x = [1, 2, 3, 4, 5, 6]
#y = [0, 0, 0, 0, 0, 0]
#z = y
#s = [.5, .6, .7, .8, .9, 1]

x = [15e3]
y = [0]
z = [0]
s = [.2]

mlab.points3d(x, y, z, s, scale_factor=1 ) #scale_mode='none'

# quiver = [x,y,z, u,v,w]
mlab.quiver3d(0, 0, 0,  r_km+10e3, 0, 0, scale_factor=1, color=(1,0,0), mask_points=5)
mlab.quiver3d(0, 0, 0,  0, r_km+10e3, 0, scale_factor=1, color=(0,1,0), mask_points=5)
mlab.quiver3d(0, 0, 0,  0, 0, r_km+10e3, scale_factor=1, color=(0,0,1), mask_points=5)


# compute magnetic field
geo = geomag.geomag.GeoMag()

lat_deg = 80.0 #
lon_deg = 190.0
alt_km = 500.0

lat = lat_deg*pi/180. #rad 
lon = lon_deg*pi/180. #rad
alt = alt_km  #km, use km as default

#deg, km
mag = geo.GeoMag(lat_deg,lon_deg,alt)

bx = mag.bx #N
by = mag.by #E
bz = mag.bz #D

print "lat:{} deg".format(mag.lat)
print "lon:{} deg".format(mag.lon)
print "alt:{} km".format(mag.alt)
#print "bh:{} nT".format(mag.bh)
print "N (bx):{} nT".format(bx)
print "E (by):{} nT".format(by)
print "D (bz):{} nT".format(bz)
#print "bz:{} nT".format(mag.bz)

total_b = np.sqrt(bx**2 + by**2 + bz**2)
print "total:{} nT".format(total_b)


# python
xx = cos(lat*pi/180.) * cos(lon*pi/180.) * (r_km + alt_km)
yy = cos(lat*pi/180.) * sin(lon*pi/180.) * (r_km + alt_km)
zz = sin(lat*pi/180.) * (r_km + alt_km) # z is 'up'

#xx = cos(phi)*cos(theta)*r
#print xx


# convert form NED to ECEF
n = bx/total_b
e = by/total_b
d = bz/total_b
##a = earth_shape; % call earth_shape to get earth data 
#a = 6378137
#b = 6356752.3142
#lat=n*pi/180
#lon=e*pi/180
#f = (a-b)/a
#e = sqrt(f*(2.-f))
#N = a /(sqrt(1.-e**2*(sin(lat))**2))
##lat = ned.pos(1)/b + ned.geo_ref(1)
#coslat = cos(lat)
#sinlat = sin(lat)
#coslon = cos(lon)
#sinlon = sin(lon)
#bx = (N + d)*coslat*coslon
#by = (N + d)*coslat*sinlon # compute the current longitude 
#bz = (N*(1 - e**2)+ d)*sinlat
##r0 = r0 + ned.geo_ref(3)*coslat
##ecef_pos.x = x # assign positions 
##ecef_pos.y = y
##ecef_pos.z = z

# Constants defined by the World Geodetic System 1984 (WGS84)
a = 6378.137
b = 6356.7523142
esq = 6.69437999014 * 0.001
e1sq = 6.73949674228 * 0.001
f = 1 / 298.257223563

## this function is ok!
def geodetic2ecef(lat, lon, alt):
    # lat, lon in radians
    # alt in km
    """Convert geodetic coordinates to ECEF."""
    #lat, lon = radians(lat), radians(lon)
    xi = sqrt(1 - esq * sin(lat))
    x = (a / xi + alt) * cos(lat) * cos(lon)
    y = (a / xi + alt) * cos(lat) * sin(lon)
    z = (a / xi * (1 - esq) + alt) * sin(lat)
    
    # x,y,z in km!!! change to meters!!! esq
    return x, y, z
    
def enu2ecef(lat, lon, alt, n, e, d):
    # lat, lon in radians
    # alt in km
    """NED (north/east/down) to ECEF coordinate system conversion."""
    x, y, z = e, n, -d
    #lat, lon = radians(lat), radians(lon)
    X, Y, Z = geodetic2ecef(lat, lon, alt)
    mx = mat('[%f %f %f; %f %f %f; %f %f %f]' %
        (-sin(lon), -sin(lat) * cos(lon), cos(lat) * cos(lon), 
          cos(lon), -sin(lat) * sin(lon), cos(lat) * sin(lon), 
          0, cos(lat), sin(lat)))
    enu = mat('[%f; %f; %f]' % (x, y, z))
    geo = mat('[%f; %f; %f]' % (X, Y, Z))
    res = mx * enu + geo
    return float(res[0]), float(res[1]), float(res[2])

def ned2ecef(lat, lon, alt, n, e, d):
    # lat, lon in radians
    #x, y, z = n, e, d
    #lat, lon = radians(lat), radians(lon)
    X, Y, Z = geodetic2ecef(lat, lon, alt)
    
    slat = sin(lat) #theta
    clat = cos(lat) #theta
    slon = sin(lon) #phi
    clon = cos(lon) #phi
    
    mx = mat('[%f %f %f; %f %f %f; %f %f %f]' %
        (-slat*clon, -slon, -clat*clon,
         -slat*slon, clon, -clat*slon,
         clat, 0, -slat)
        )
    ned = mat('[%f; %f; %f]' % (n, e, d))
    geo = mat('[%f; %f; %f]' % (X, Y, Z))
    
    res = mx * ned
    return float(res[0]), float(res[1]), float(res[2])




#def plotMagVector(lat,lon,alt):
    
bx_ecef, by_ecef, bz_ecef = ned2ecef(lat,lon,alt,n,e,d)
xx_ecef, yy_ecef, zz_ecef = geodetic2ecef(lat,lon,alt)
    #print bx_ecef,by_ecef,bz_ecef
    #print xx_ecef,yy_ecef,zz_ecef
    #print xx,yy,zz
mlab.quiver3d(  xx_ecef, 
                yy_ecef, 
                zz_ecef, 
                bx_ecef, 
                by_ecef, 
                bz_ecef, 
                scale_factor=3000, color=(0,0,1))
    
#plotMagVector(lat,lon,alt)

#print ned2ecef(0,90*pi/180,0,7,0,0)
#mlab.view(63.4, 73.8, 4, [-0.05, 0, 0])
#mlab.show()


