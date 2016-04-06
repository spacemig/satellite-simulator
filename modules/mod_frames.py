from numpy import linalg, array, cross
from math import pow, sin, cos, sqrt
from scipy import mat, arctan, arctan2 #pi, cos, sin,
from mod_constants import wgs_constants, vector_zero 

## this function is ok!
def ecef_from_geodetic(lat, lon, alt):
    # load wgs constants
    wgs = wgs_constants()
    # lat, lon in radians
    # alt in km
    """Convert geodetic coordinates to ECEF."""
    #lat, lon = radians(lat), radians(lon)
    xi = sqrt(1 - wgs.esq * sin(lat))
    x = (wgs.a / xi + alt) * cos(lat) * cos(lon)
    y = (wgs.a / xi + alt) * cos(lat) * sin(lon)
    z = (wgs.a / xi * (1 - wgs.esq) + alt) * sin(lat)
    
    # x,y,z in km!!! change to meters!!! esq
    return x, y, z
    
def ecef_from_enu(lat, lon, alt, n, e, d):
    # lat, lon in radians
    # alt in km
    """NED (north/east/down) to ECEF coordinate system conversion."""
    x, y, z = e, n, -d
    #lat, lon = radians(lat), radians(lon)
    X, Y, Z = ecef_from_geodetic(lat, lon, alt)
    mx = mat('[%f %f %f; %f %f %f; %f %f %f]' %
        (-sin(lon), -sin(lat) * cos(lon), cos(lat) * cos(lon), 
          cos(lon), -sin(lat) * sin(lon), cos(lat) * sin(lon), 
          0, cos(lat), sin(lat)))
    enu = mat('[%f; %f; %f]' % (x, y, z))
    geo = mat('[%f; %f; %f]' % (X, Y, Z))
    res = mx * enu + geo
    return float(res[0]), float(res[1]), float(res[2])

def ecef_from_ned(lat, lon, alt, n, e, d):
    # lat, lon in radians
    #x, y, z = n, e, d
    #lat, lon = radians(lat), radians(lon)
    #X, Y, Z = geodetic2ecef(lat, lon, alt)
    
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
    #geo = mat('[%f; %f; %f]' % (X, Y, Z))
    
    res = mx * ned
    return float(res[0]), float(res[1]), float(res[2])
    
def cbrt(x):
    if x >= 0: 
        return pow(x, 1.0/3.0)
    else:
        return -pow(abs(x), 1.0/3.0)
           
def geodetic_from_ecef(x, y, z):
    #http://code.google.com/p/pysatel/source/browse/trunk/coord.py?r=22
    """Convert ECEF coordinates to geodetic.
    J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates \
    to geodetic coordinates," IEEE Transactions on Aerospace and \
    Electronic Systems, vol. 30, pp. 957-961, 1994."""

    # load wgs constants
    wgs = wgs_constants()
    a = wgs.a
    b = wgs.b
    esq = wgs.esq
    e1sq = wgs.e1sq
    
    r = sqrt(x * x + y * y)
    Esq = a * a - b * b
    F = 54 * b * b * z * z
    G = r * r + (1 - esq) * z * z - esq * Esq
    C = (esq * esq * F * r * r) / (pow(G, 3))
    S = cbrt(1 + C + sqrt(C * C + 2 * C))
    P = F / (3 * pow((S + 1 / S + 1), 2) * G * G)
    Q = sqrt(1 + 2 * esq * esq * P)
    r_0 =  -(P * esq * r) / (1 + Q) + sqrt(0.5 * a * a*(1 + 1.0 / Q) - \
        P * (1 - esq) * z * z / (Q * (1 + Q)) - 0.5 * P * r * r)
    #U = sqrt(pow((r - esq * r_0), 2) + z * z)
    V = sqrt(pow((r - esq * r_0), 2) + (1 - esq) * z * z)
    Z_0 = b * b * z / (a * V)
    #h = U * (1 - b * b / (a * V))
    lat = arctan((z + e1sq * Z_0) / r)
    lon = arctan2(y, x)
    return lat, lon
    #return degrees(lat), degrees(lon)

#geodetic_from_ecef(1e6,2e6,3e6)


def body_frame(position,velocity):
    #compute the orbital frame from the position and velocity
    
    position0 = position
    velocity0 = velocity
    
    body_x = array(velocity0)
    body_x = body_x/linalg.norm(body_x)
    body_x = body_x.tolist()
    
    #vector_zero = [0,0,0]
    
    body_z = vector_zero - array(position0)
    
    # [i+j for i,j in zip(L, M)]
    #body_z = [i-j for i,j in zip(vector_zero, position0)]
    
    # norm vector
    body_z = body_z/linalg.norm(body_z)
    body_z = body_z.tolist()
    
    # y axis = cross(z,x)
    
    body_y = cross(body_z,body_x)
    body_y = body_y.tolist()
    
    return body_x,body_y,body_z