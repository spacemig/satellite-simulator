from numpy import *

## gravitational parameter
# miu = G*(m_earth + m_satellite) = G*(m_earth)
G = 6.673e-20 #+/- 0.001e-20 km^3/(kg.s^2)

## Mass of the Earth
m_earth = 5.9733320e24 # [kg]

## Earth mean radius 
r_earth = 6378.1363 #km

## miu = G*m_earth
miu = 398600.4415 #km^3/s^2

## satellite altitude above the surface of the Earth
# for circular orbits
h = 590.0 #[km]

## satellite altitude
r_satellite= r_earth + h

# angular rate [rad/s]
omega = sqrt(miu/r_satellite**3)

# gravity acceleration
g = miu/(r_satellite)**2*1e3 # convert to m/s^2 (not km/s^2)

# motor specifics
Isp = 320 #[s]

c = Isp*g

class wgs_constants:
    # Constants defined by the World Geodetic System 1984 (WGS84)
    a = 6378.137
    b = 6356.7523142
    esq = 6.69437999014 * 0.001
    e1sq = 6.73949674228 * 0.001
    f = 1 / 298.257223563

vector_zero = array([0,0,0])