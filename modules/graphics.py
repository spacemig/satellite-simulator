#graphics support library

from mayavi import mlab
import numpy as np
import orbit
import frames
from tvtk.tools import visual
from math import pow, degrees, radians, pi, sin, cos, sqrt

vector_zero = orbit.vector_zero

def drawVector(origin,vector,scale=1,color=(1,0,0)):
    mlab.quiver3d(  
                origin[0], 
                origin[1], 
                origin[2], 
                vector[0], 
                vector[1], 
                vector[2], 
                scale_factor=scale, color=color)
      
def drawContinentsSphere(r_earth):          
    ###############################################################################
    # Display continents outline, using the VTK Builtin surface 'Earth'
    from mayavi.sources.builtin_surface import BuiltinSurface
    continents_src = BuiltinSurface(source='earth', name='Continents')
    # The on_ratio of the Earth source controls the level of detail of the
    # continents outline.
    #r_earth = 6371.0 # km 
    continents_src.data_source.on_ratio = 2
    continents_src.data_source.radius = r_earth
    continents = mlab.pipeline.surface(continents_src, color=(0, 0, 0))


def drawGlobe(r_earth):
    ###############################################################################
    # Display a semi-transparent sphere, for the surface of the Earth
    
    # We use a sphere Glyph, throught the points3d mlab function, rather than
    # building the mesh ourselves, because it gives a better transparent
    # rendering.
    ocean_blue = (0.4, 0.5, 1.0)
    sphere = mlab.points3d(0, 0, 0, scale_mode='none',
                                    scale_factor=r_earth*2,
                                    color=ocean_blue,
                                    resolution=50,
                                    opacity=0.5,
                                    name='Earth')
    
    # These parameters, as well as the color, where tweaked through the GUI,
    # with the record mode to produce lines of code usable in a script.
    sphere.actor.property.specular = 0.45
    sphere.actor.property.specular_power = 5
    # Backface culling is necessary for more a beautiful transparent
    # rendering.
    sphere.actor.property.backface_culling = True
    
def drawEquator(r_earth):
    ###############################################################################
    # Plot the equator and the tropiques
    theta = np.linspace(0, 2*np.pi, 100)
    for angle in (-np.pi/6, 0, np.pi/6):
        x = r_earth*np.cos(theta)*np.cos(angle)
        y = r_earth*np.sin(theta)*np.cos(angle)
        z = r_earth*np.ones_like(theta)*np.sin(angle)
    
        mlab.plot3d(x, y, z, color=(1, 1, 1),
                            opacity=0.1, tube_radius=None)
    
###############################################################################                        
def drawReferenceFrameEarth(r_earth): 
    """ draw x, y z axes
    """
    # using quiver = [x,y,z, u,v,w]
    mlab.quiver3d(0, 0, 0,  r_earth+10e3, 0, 0, scale_factor=1, color=(1,0,0), mask_points=5)
    mlab.quiver3d(0, 0, 0,  0, r_earth+10e3, 0, scale_factor=1, color=(0,1,0), mask_points=5)
    mlab.quiver3d(0, 0, 0,  0, 0, r_earth+10e3, scale_factor=1, color=(0,0,1), mask_points=5)

def drawReferenceFrameFromVectors(pos,v1,v2,v3,size):
    #print pos
    # check if Visual.Mframe can draw reference frames

    """ This function draws 3 orthogonal vectors to represent a reference frame 
    with X, Y, Z axes.
    
    :param pos: the position array (x,y,z) in meters
    :param v1: vector 1
    :param v2: vector 2
    :param v3: vector 3
    :returns: mlab quiver objects
    
    """

    #pos = pos.tolist()
    
    # quiver usage = [x,y,z, u,v,w]
    a = mlab.quiver3d(pos[0], pos[1], pos[2],  v1[0], v1[1], v1[2], scale_factor=size, color=(1,0,0), mask_points=5)
    b = mlab.quiver3d(pos[0], pos[1], pos[2],  v2[0], v2[1], v2[2], scale_factor=size, color=(0,1,0), mask_points=5)
    c = mlab.quiver3d(pos[0], pos[1], pos[2],  v3[0], v3[1], v3[2], scale_factor=size, color=(0,0,1), mask_points=5)
    return [a,b,c]
    
def drawReferenceFrameFromDCM(pos,dcm,size):
    # check Visual.Mframe
    ###############################################################################
    # draw x, y z axes
    # quiver = [x,y,z, u,v,w]
    
    v1 = np.dot(dcm,np.array([1,0,0]))
    v2 = np.dot(dcm,np.array([0,1,0]))
    v3 = np.dot(dcm,np.array([0,0,1]))
    
    pos = pos.tolist()
    a = mlab.quiver3d(pos[0], pos[1], pos[2],  v1[0], v1[1], v1[2], scale_factor=size, color=(1,0,0), mask_points=5)
    b = mlab.quiver3d(pos[0], pos[1], pos[2],  v2[0], v2[1], v2[2], scale_factor=size, color=(0,1,0), mask_points=5)
    c = mlab.quiver3d(pos[0], pos[1], pos[2],  v3[0], v3[1], v3[2], scale_factor=size, color=(0,0,1), mask_points=5)
    return [a,b,c]

def drawMagVector(lat,lon,alt):
    # lat and lon in deg
    # plot and compute magnetic field
    #geo = geomag.geomag.GeoMag()
    #lat = lat_deg*np.pi/180. #rad 
    #lon = lon_deg*np.pi/180. #rad
    #alt = alt_km  #km, use km as default

  
    bx_ecef, by_ecef, bz_ecef, b_total = orbit.mag_EFEC(lat,lon,alt)
    xx_ecef, yy_ecef, zz_ecef = orbit.geodetic2ecef(lat,lon,alt)
    #print bx_ecef,by_ecef,bz_ecef
    #print xx_ecef,yy_ecef,zz_ecef
    
    b_ecef_unit = orbit.normalize(np.array([bx_ecef, by_ecef, bz_ecef]))

    mlab.quiver3d(  xx_ecef, 
                yy_ecef, 
                zz_ecef, 
                float(b_ecef_unit[0]), 
                float(b_ecef_unit[1]), 
                float(b_ecef_unit[2]), 
                scale_factor=2000, color=(0,1,1))
                
def drawSatellitePosition(dateTime,tle):
    # dateTime = year, month, day, hour, minute, second
    position, velocity = orbit.state_satellite(dateTime,tle)
    x = position[0]
    y = position[1]
    z = position[2]
    
    #x, y, z = 10e3,0,0
    # draw satellite position
    s = [200] #scale
    mlab.points3d(x, y, z, s, scale_factor=1, color=(1,1,0) ) #scale_mode='none'

def drawSatellite(position,euler):

    # euler angles are in radians
    
    psi, theta, phi = euler
    
    sat = visual.cylinder()
    sat.length = 300
    sat.radius = 300
    sat.x = position[0]
    sat.y = position[1]
    sat.z = position[2]
    
    sat.rotate(degrees(psi)   , [0,0,1] ,position)
    sat.rotate(degrees(theta) , [0,1,0] ,position)
    sat.rotate(degrees(phi)   , [1,0,0] ,position)

    return sat

def drawSatelliteState(satobj,position,euler):

    # euler angles are in radians
    
    psi, theta, phi = euler
    
    #sat = visual.cylinder()
    
    sat = satobj
    sat.length = 300
    sat.radius = 300
    sat.x = position[0]
    sat.y = position[1]
    sat.z = position[2]
    
    sat.rotate(degrees(psi)   , [0,0,1] ,position)
    sat.rotate(degrees(theta) , [0,1,0] ,position)
    sat.rotate(degrees(phi)   , [1,0,0] ,position)
    
    return 1

def drawOrbitalFrame(position,velocity):
    # Inputs:
    # - position, satellite position
    # - velocity, satellite velocity
    ###############################################################################
    # plot orbital frame
    #mlab.quiver3d(  position0[0], 
    #                position0[1], 
    #                position0[2], 
    #                velocity[0], 
    #                velocity[1], 
    #                velocity[2], 
    #                scale_factor=500, color=(0,0,1))
    
    position0 = position
    velocity0 = velocity
    
    #body_x, body_y, body_z = orbit.bodyFrame(position,velocity)
    body_x, body_y, body_z = frames.body_frame(position,velocity)
    
    # X axis
    drawVector(position0,body_x,3000,(1,0,0))
    #mlab.quiver3d(  position0[0], position0[1], position0[2],  
    #                body_x[0], body_x[1], body_x[2], scale_factor=3000, color=(1,0,0), mask_points=5)
    

    # Y axis
    #float(body_z[0])
    drawVector(position0,body_y,3000,(0,1,0))
    #mlab.quiver3d(  position0[0], position0[1], position0[2],  
    #                body_z[0], body_z[1], body_z[2], scale_factor=3000, color=(0,1,0), mask_points=5)
    
    # Z axis
    drawVector(position0,body_z,3000,(0,0,1))
    #mlab.quiver3d(  position0[0], position0[1], position0[2],  
    #                body_y[0], body_y[1], body_y[2], scale_factor=3000, color=(0,0,1), mask_points=5)


def drawOrbit(dateTime_,nk,tle):
    # dateTime = [year, month, day, hour, minute, second]
    ###############################################################################
    # Compute Orbit Points and Plot them
    position = []
    # change to just dateTime
    position0, velocity0 = orbit.state_satellite(dateTime_,tle)
    
    # I have to convert back to individual values because there is a strange 
    # error when I do dateTime_[4] = dateTime_[4]+k ... this keeps incrementing
    # the original dateTime_
    year, month, day, hour, minute, second = dateTime_
    #dateTime = dateTime_
    
    for k in range(0,nk):
        # step for each minute
        # collect all positions for satellite
        
        # increment the minutes
        #dateTime[4] = dateTime_[4]+k # [year, month, day, hour, minute+k, second]
        position_, velocity_ = orbit.state_satellite([year, month, day, hour, minute+k, second],tle)
        #position_, velocity_ = orbit.state_satellite(dateTime)
        
        position.append(position_)
        
        #position = np.append([position], [position0], axis=0)
        
        #return position
        #x0 = position0[0]
        #y0 = position0[1]
        #z0 = position0[2]
    
    # convert list to array
    position = np.array(position)
    
    x = position[:,0]
    y = position[:,1]
    z = position[:,2]
    
    
    #position0, velocity0 = state_satellite(year, month, day, hour, minute, second)
    #x0 = position0[0]
    #y0 = position0[1]
    #z0 = position0[2]
    #
    #position1, velocity1 = state_satellite(year, month, day, hour, minute+10, second)
    #x1 = position1[0]
    #y1 = position1[1]
    #z1 = position1[2]
    
    # draw orbit
    #x = np.array([x0, x1])
    #y = np.array([y0, y1])
    #z = np.array([z0, z1])
    #mlab.plot3d(x,y,z,  tube_radius=50, color=(0,0,1)) #colormap='Spectral'
    
    # return the mlab object in case of wanting to manipulate it
    return mlab.plot3d(x,y,z,  tube_radius=10, color=(0,0,1)), position