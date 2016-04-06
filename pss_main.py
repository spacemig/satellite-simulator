#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Miguel Nunes"
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Miguel Nunes"
__email__ = "spacemig@gmail.com"
__status__ = "Production"

# Python Space Simulator v 0.0.1

# make sure you have the following modules installed
# - geomag
# for help on installing these 

import os
import sys
# system modules
from time import time, sleep
tic = time()

#print __file__
#print("Path at terminal when executing this file")
#print(os.getcwd() + "\n")
dirname, filename = os.path.split(os.path.abspath(__file__))
#print dirname

# get the current script file path
this_folder      = os.path.dirname(os.path.abspath(__file__))+'/'
modules_folder   = this_folder+'modules/'
geomag_folder    = this_folder+'thirdparty/geomag-0.9/'
sgp4_folder      = this_folder+'thirdparty/sgp4/'
pyorbital_folder = this_folder+'thirdparty/pyorbital/'
tle_folder       = this_folder+'thirdparty/tle/'

# add the folders to the system path so we don't have to install the folders
# althought that is recommended
if modules_folder not in sys.path:
    sys.path.append(modules_folder)
    
if geomag_folder not in sys.path:
    sys.path.append(geomag_folder)

if sgp4_folder not in sys.path:
    sys.path.append(sgp4_folder)

if pyorbital_folder not in sys.path:
    sys.path.append(pyorbital_folder)

if tle_folder not in sys.path:
    sys.path.append(tle_folder)

import mod_orbit_dynamics as osl #orbital dynamics support library
import mod_graphics as gsl #graphics support library
import mod_attitude as attitude

# thirparty modules
from pyorbital import tlefile
from pyorbital.orbital import Orbital
from datetime import datetime

# math modules
import numpy as np
from numpy import dot, array, shape
from numpy.linalg import norm
from math import pow, degrees, radians, pi, sin, cos, sqrt
from scipy import mat,  arctan, arctan2, integrate #pi, cos, sin,

# graphics modules
from mayavi import mlab
import matplotlib.pyplot as plt
from pylab import ion #, plots
from tvtk.tools import visual



# $python setup.py install for first time
#import geomag 

# -----------------------------------------------------------------------------
# Load the ISS TLE
#tle = tlefile.read('ISS (ZARYA)') # reads TLE from the celestrack website
tle = tlefile.read('ISS (ZARYA)',tle_folder+'stations.txt')
#print tle.epoch

orb = Orbital("ISS (ZARYA)",tle_folder+'stations.txt')
#now = datetime.utcnow() # for real time
now = datetime(2013,05,31,22,0,0)
pos0 = orb.get_position(now,normalize=False)
lon,lat,alt = orb.get_lonlatalt(now)
print "\nSatellite Lat Lon Alt: ",lat,lon,alt,"\n"
print "\nSatellite Position (normalized):", pos0[0], "\n"
vector_zero = osl.vector_zero

###############################################################################
r_earth = 6371.0 # km 
r = r_earth

# start epoch
year = 2013
month = 11
day = 1
hour = 0
minute = 30
second = 0

year = now.year
month = now.month
day = now.day
hour = now.hour
minute = now.minute
second = now.second


dateTime = [year, month, day, hour, minute, second]

# mean motion (revs/day)
mean_motion = 16.10627010
orbit_period = (60*60*24)/mean_motion

###############################################################################
# MayaVi figure

# Create a figure
figMaya = mlab.figure(1, 
    bgcolor=(0.48, 0.48, 0.48), 
    fgcolor=(0, 0, 0),
    size=(400, 400))

# remove all actors in the figure
n_actors = figMaya.scene.renderer.actors.number_of_items
for k in range(n_actors):
    #print k
    figMaya.scene.renderer.remove_actor(figMaya.scene.renderer.actors[0])


#mlab.clf()

# Tell visual to use this as the viewer.
visual.set_viewer(figMaya)

gsl.drawContinentsSphere(r_earth)
gsl.drawGlobe(r_earth)
gsl.drawEquator(r_earth)
gsl.drawReferenceFrameEarth(r_earth)


#lat_deg = 80.0 #
#lon_deg = 190.0
#alt_km = 500.0

#plotMagVector(-50,50,500)

# plot a range of magnetic vectors
#for k in range(-90,90,10):
#    gsl.drawMagVector(k,0,500)

#gsl.drawSatellitePosition(dateTime,tle)
#b1 = visual.box()
#b1.size = (1000,1000,1000)

position_sat, velocity_sat = osl.state_satellite(dateTime,tle)

print "\nSatellite Position: ", position_sat, "\n"
#position_sat = [-6479.7,-1762,976]
#c1 = visual.cylinder()
#c1.length = 300
#c1.radius = 300
#c1.x = position_sat[0]
#c1.y = position_sat[1]
#c1.z = position_sat[2]

#c2 = visual.cylinder()
#c2.length = 300
#c2.radius = 300
#c2.x = 7000
#c2.y = 0
#c2.z = 0

#sat = gsl.drawSatellite(position_sat,[0,0,0])

orbit,orbit_position = gsl.drawOrbit(dateTime,15,tle)
pos = array([position_sat[0],position_sat[1],position_sat[2]])
#fr = gsl.drawReferenceFrame(pos,2000)
ofr = gsl.drawOrbitalFrame(position_sat,velocity_sat)

# testing rotations
a = array([7000,0,0])
#gsl.drawVector(vector_zero,a)

# yaw-psi, pitch-theta, roll-phi
yaw,pitch,roll = radians(45),radians(90),radians(0)
psi = yaw
theta = pitch
phi = roll 

#visual._create_rotation_matrix angle in degrees, 
# these are active transformations 
Rx = visual._create_rotation_matrix([1,0,0],degrees(phi))
Ry = visual._create_rotation_matrix([0,1,0],degrees(theta))
Rz = visual._create_rotation_matrix([0,0,1],degrees(psi))

Rx = attitude.rotX3d_passive(phi)
Ry = attitude.rotY3d_passive(theta)
Rz = attitude.rotZ3d_passive(psi)
# aerospace standard
# first rotate around z, then y, then x
# Rx*Ry*Rz
Rtotal = dot(Rz,dot(Ry,Rx))

#b1 = dot(Rx,dot(Ry,dot(Rz,a)))
b = dot(Rtotal.T,a)
gsl.drawVector(vector_zero,b)
#gsl.drawVector(vector_zero,b2)


#q1 = attitude.euler2quaternion(radians(r),radians(p),radians(y))
#euler = attitude.quaternion2euler(q1)*180/pi

#print euler
#print q1

#DCM = attitude.quaternion2dcm(q1)
#print DCM
#print Rtotal

# this works
#q1 = attitude.euler2quaternion_test(array([radians(r),radians(p),radians(y)]))
#euler1 = attitude.quaternion2euler_test(array(q1))*180/pi
euler = array([psi,theta,phi]) # use, 
q1 = attitude.quaternion_from_euler(euler)
euler1 = attitude.euler_from_quaternion(array(q1))*180/pi

print q1
print euler1

DCM = attitude.dcm_from_quaternion(q1)
print Rtotal,"\n"
print DCM,"\n"
print attitude.dcm_from_euler(euler),"\n"
a  = array([1,0,0])
a_ = dot(DCM,a)

b  = array([1,0,0])
b_ = dot(Rtotal,b)

#########################################
# Initial conditions
#########################################

time_data = range(0,1000,1)

# define initial state vector
roll  = 50*pi/180 # roll in deg
pitch = 20*pi/180 # pithc in deg
yaw   = 0*pi/180 # yaw in deg

#q_0 = tr.quaternion_from_euler(roll,pitch,yaw)
euler_0 = array([roll,pitch,yaw])
q_0 = attitude.quaternion_from_euler(euler_0)
#q_0 = array([0,0,0,1])
omega_x0 = 2*pi/180
omega_y0 = 1*pi/180
omega_z0 = 0*pi/180

omega_0 = array([omega_x0, omega_y0, omega_z0])
#x_0 = [0.0, 0.0, 0.2, 1.0, 0.0, 0.0, 0.0]  #intial state vector: omega_x0,omega_y0,omega_z0,q0,q1,q2,q3  
#x_0 = omega_0 + q_0.tolist()
x_0 = np.hstack((omega_0,q_0))

torque_0 = array([0.0, 0.0, 0.0]) # initial torque
#
##########################################
##a = array ([1,2,3])
##B = array( [ (1,2,3), (4,5,6), (7,8,9) ] )
#
#Inertia = array( [(2.4388, 0.0573, -0.0396), 
#            (0.0573, 2.4155, -0.0260), 
#            (-0.0396, -0.0260, 2.8271) ] )
#
#using principal moments of inertia
I_x = 2.4388
I_y = 2.4155
I_z = 2.8271

Inertia = array( [(I_x   , 0     , 0     ), 
                  (0     , I_y   , 0     ), 
                  (0     , 0     , I_z) ] )

#########################################
# SIM
#########################################

x_sim = integrate.odeint(attitude.attitude_dynamics, x_0, time_data, args=(torque_0,Inertia))

## check the quaternion module
##q_mod = sqrt(x_sim[:,3]**2+x_sim[:,4]**2+x_sim[:,5]**2+x_sim[:,6]**2)
##print x_sim
#
## assign results 
#omega_sim = x_sim[:,0:3]
#q_sim = x_sim[:,3:7]
#
## convert all quaternions to euler angles
#rpy = np.zeros((q_sim.shape[0],3))
#for i in range(0,q_sim.shape[0]):
#    rpy[i,:] = attitude.quaternion2euler_test(q_sim[i,:])

##########################################
## INIT FIGURE
##########################################
#
##turn interactive mode on
#ion()
#
## to get size: fig.get_size_inches()
#fig = plt.figure(1) #figsize=(6.5, 9.5))
#
##clear the figure
#fig.clf()
##fig.canvas.manager.window.move(900,0) # in pixels from top-left corner
#
##ax = fig.add_subplot(111)
#
#
#
##########################################
## PLOTS
##########################################

## plot measurements and compare with estimated fiter measurements
#plt.subplot(311)
#plt.plot(time_data, omega_sim[:,0]) # t, omega_meas*180/pi,
##ylim(omega_ekf.min()*180/pi*1.1-1, omega_ekf.max()*180/pi*1.1+1)
##legend(['$\omega_x$','$\omega_y$','$\omega_z$'])
##ylabel('$\Omega$ [deg/sec]')
#plt.subplot(312)
#plt.plot(time_data, omega_sim[:,1]) 
#
#plt.subplot(313)
#plt.plot(time_data, omega_sim[:,2]) 
#
#
#fig = plt.figure(2)
#fig.clf()
#plt.subplot(411)
#plt.plot(time_data, q_sim[:,0]) 
#
#plt.subplot(412)
#plt.plot(time_data, q_sim[:,1]) 
#
#plt.subplot(413)
#plt.plot(time_data, q_sim[:,2]) 
#
#plt.subplot(414)
#plt.plot(time_data, q_sim[:,3]) 
#
####
#fig = plt.figure(3)
#fig.clf()
#plt.subplot(311)
#plt.plot(time_data, rpy[:,0]*180/pi) 
#plt.title('Euler Angles')
#
#plt.subplot(312)
#plt.plot(time_data, rpy[:,1]*180/pi) 
#
#plt.subplot(313)
#plt.plot(time_data, rpy[:,2]*180/pi) 

'''
###############################################################################
# plot magnetic vector in satellite position
#

lat, lon = osl.ecef2geodetic(position0[0],position0[1],position0[2])
alt = float(norm(np.array(position0)) - r_earth)
#print alt
gsl.drawMagVector(lat,lon,alt)
#gsl.drawOrbitalFrame(position0,velocity0)

##############################################################################
# testing draw mag vector 
#lat,lon,alt = 0,0,0
#bx_ecef, by_ecef, bz_ecef, b_total = osl.mag_EFEC(lat,lon,alt)
#xx_ecef, yy_ecef, zz_ecef = osl.geodetic2ecef(lat,lon,alt)
#mlab.quiver3d( 7000,0,0, 1,0,0, scale_factor=7000, color=(0,1,1))
#mlab.quiver3d( 7000,0,0, 1,1,0, scale_factor=7000, color=(0,1,1))

#print lat,lon,alt
#bx_ecef_unit = osl.normalize(array([bx_ecef, by_ecef, bz_ecef]))

#print(bx_ecef_unit)

#print bx_ecef, by_ecef, bz_ecef, norm(array([bx_ecef, by_ecef, bz_ecef]))*1e9
#print xx_ecef, yy_ecef, zz_ecef
#mlab.quiver3d( 7000,0,0, float(bx_ecef_unit[0]), float(bx_ecef_unit[1]), float(bx_ecef_unit[2]), scale_factor=7000, color=(0,1,1))


#print ned2ecef(0,90*pi/180,0,7,0,0)
#mlab.view(63.4, 73.8, 4, [-0.05, 0, 0])
#mlab.show()

#mlab.view(0,90,5000,[position0[0], position0[1], position0[2]])

# compute the angle between the body x vector and the magnetic vector
# this will help to determine the position of the torque rods
#body_x = osl.normalize(np.array(velocity0))
body_x, body_y, body_z = osl.bodyFrame(position0,velocity0)



#body_x = np.array(body_x)
bx_ecef, by_ecef, bz_ecef, b_total = osl.mag_EFEC(lat,lon,alt)
#print bx_ecef,by_ecef,bz_ecef
mag_ = np.array([bx_ecef, by_ecef, bz_ecef])

# angle between two vectors
#angle = np.arccos(np.dot(body_x,mag_)/(norm(body_x)*norm(mag_)))*180/pi
#print angle


# DCM for <body_x,X>
#Frame A, EFEC
X = [1,0,0]
Y = [0,1,0]
Z = [0,0,1]

#Frame B
#body_x = [1,1,0]
#body_y = [-1,1,0]
#body_z = [0,0,1]

dcos = osl.dcos
def dcm(ax,ay,az, bx,by,bz):
    return np.array([[ dcos(bx,ax), dcos(bx,ay), dcos(bx,az)],
                    [  dcos(by,ax), dcos(by,ay), dcos(by,az)],
                    [  dcos(bz,ax), dcos(bz,ay), dcos(bz,az)]]) 
      
                            
#DCM = np.array([[ dcos(body_x,X), dcos(body_x,Y), dcos(body_x,Z)],
#                 [dcos(body_y,X), dcos(body_y,Y), dcos(body_y,Z)],
#                 [dcos(body_z,X), dcos(body_z,Y), dcos(body_z,Z)]]) 
DCM = dcm(X,Y,Z, body_x,body_y,body_z)

#convert a vector in the EFEC frame to the body frame (BF)   
# this computes the EFEC coordinates of the X vector rotated to the BF, from EFEC orientation
#newX = dot(DCM.T,X)
# what I want is the BF coordintates of some vector in the EFEC 
# this computes the BF coordinates of the X vector resident in the EFEC
#newXX = dot(DCM,X)
#gsl.drawVector(vector_zero,array(X),1,(1,1,0))
#gsl.drawVector(array(position0),newX,1,(1,1,0))   
#print newX 
#print newXX

#convert the Magnetic Vector in the EFEC to the BF
#newMag = dot(DCM,mag_)
#gsl.drawVector(array(position0),newMag,1000,(1,1,0))   
#print newMag

###############################################################################
## Plot the deviation of the magnetic vector with the body x axis (velocity vector)
position_data = []
torque_data = []
#position0, velocity0 = osl.state_satellite(year, month, day, hour, minute, second)

time_data = range(0,10,2)

# magnetic moment, maximum for trods
# http://www.vectronic-aerospace.com/space.php?p=torquer
magnetic_moment_max = 25 #A.m^2

for k in time_data:
    # step for each minute
    # collect all positions for satellite
    position_, velocity_ = osl.state_satellite(year, month, day, hour, minute+k, second)
    
    position_data.append(position_)
    
    lat, lon = osl.ecef2geodetic(position_[0],position_[1],position_[2])
    alt = float(norm(np.array(position_)) - r_earth)
    bx_ecef, by_ecef, bz_ecef, b_total = osl.mag_EFEC(lat,lon,alt)
    mag_ = np.array([bx_ecef, by_ecef, bz_ecef])
    body_x, body_y, body_z = osl.bodyFrame(position_,velocity_)
    
    DCM = dcm(X,Y,Z, body_x,body_y,body_z)
    
    #convert the Magnetic Vector coodinates in the EFEC to the BF coordinates
    mag_bf = dot(DCM,mag_)
    #print mag_bf
    
    # now compute the angle between the mag_bf projected in the bx,by plane of the body axes!!!!
    #angle = arctan(mag_bf[1]/mag_bf[0])*180./pi
    #angle = np.arccos(np.dot(body_x,mag_)/(norm(body_x)*norm(mag_)))*180/pi
    #angle = np.arccos(dcos([1,0],[mag_bf[0],mag_bf[1]]))*180.0/pi
    #print angle
    
    # compute the cross product between the mag vector projection on the X,Y plane of the body
    # this gives the the total torque for yaw maneuvers from the torque rod - X 
    
    # adimensional 
    #tx_max = np.cross([1,0],[mag_bf[0],mag_bf[1]])
    #ty_max = np.cross([0,1],[mag_bf[0],mag_bf[1]])
    
    #dimensional 
    tx_max = np.cross([magnetic_moment_max,0],[mag_bf[0],mag_bf[1]])
    ty_max = np.cross([0,magnetic_moment_max],[mag_bf[0],mag_bf[1]])
    #print tx_max, ty_max
    
    torque_data.append([tx_max,ty_max])
    
    gsl.drawMagVector(lat,lon,alt)
    #gsl.drawOrbitalFrame(position_,velocity_)
    
    # REVISIT THE PREVIOUS CODE... really need position0 inside the loop?
    


#########################################
# PLOTS
#########################################    
#turn interactive mode on
ion()
fig1 = plt.figure(1)
fig1.clf() #clear the figure

plt.subplot(211)
plt.plot(array(time_data)*60./orbit_period,array(torque_data)[:,0])
plt.ylim(-0.001,0.001)
plt.ylabel('Mag Torquer X (N.m)')
plt.title('HiakaSat - Max Possible Torque from Trods for Yaw control (106 $^o$ inc.) \n (assuming trod-x is aligned with velocity vector)')

plt.subplot(212)
plt.plot(array(time_data)*60./orbit_period,array(torque_data)[:,1])
plt.ylim(-0.001,0.001)
plt.ylabel('Mag Torquer Y (N.m)')
plt.xlabel('Time (orbits)')
'''


###############################################################################
# Animation

kk = 0
jjj = 0
number_iterations = int(shape(orbit_position)[0]) -10 
end_flag = 0

def satellite_state(kk):
    #global kk
    #kk+=1
    
    #print(time_data[kk])
    #orbit_position
    q = x_sim[kk,3:7]
    euler = attitude.euler_from_quaternion(q)
    
    return orbit_position[kk,0], orbit_position[kk,1], orbit_position[kk,2],euler[0],euler[1],euler[2]

    
#@mlab.show
#@mlab.animate(delay=1000, ui=True)
#def anim():
#    """Animate the satellite """
#    #t = t+1
#    #jjj < number_iterations:
#    fgcf = mlab.gcf()
#    while 1: 
#        #jjj = jjj + 1
#        x, y, z, phi,theta,psi = satellite_state()
#        print x, y, z, degrees(phi),degrees(theta),degrees(psi)
#        #sat.x, sat.y, sat.z, phi,theta,psi = satellite_state()
#        #sat = gsl.drawSatellite([x,y,z],[phi,theta,psi])
#        gsl.drawSatelliteState(sat,[x,y,z],[phi,theta,psi])
#        
#        fgcf.scene.render()
#        
#        yield
#
### Run the animation.
#anim()

# 
for k in range(number_iterations):
    
    print k
    x, y, z, phi,theta,psi = satellite_state(k)
    sat = gsl.drawSatellite([x,y,z],[phi,theta,psi])

###############################################################################
# ELAPSED TIME
############################################################################### 
toc = time()

print "Elapsed time: {:.3f} s".format(toc-tic)