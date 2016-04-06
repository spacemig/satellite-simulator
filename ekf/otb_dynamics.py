# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 10:31:17 2013

@author: miguelnunes

references: 
    

"""


from adcs_support_lib import *
import time
start = time.time()

print "\n#######################################\n"

#########################################
# Initial conditions
#########################################

t = range(0,100,1)

# define initial state vector
roll  = 0*pi/180 # roll in deg
pitch = 0*pi/180 # pithc in deg
yaw   = 0*pi/180 # yaw in deg

#q_0 = array([0,0,0,1])
q_0 = euler2quaternion(roll,pitch,yaw)

# initial angular rates
omega_x0 = 0*pi/180
omega_y0 = 0*pi/180
omega_z0 = 10*pi/180

omega_0 = array([omega_x0, omega_y0, omega_z0])
#x_0 = [0.0, 0.0, 0.2, 1.0, 0.0, 0.0, 0.0]  #intial state vector: omega_x0,omega_y0,omega_z0,q0,q1,q2,q3  
#x_0 = omega_0.tolist() + q_0.tolist()
x_0 = array([omega_x0, omega_y0, omega_z0, q_0[0], q_0[1], q_0[2], q_0[3] ])
torque_0 = [0.0, 0.0, -0.0001] # initial torque

#########################################

#using principal moments of inertia
I_x = 2.4388
I_y = 2.4155
I_z = 2.8271

#I_x = 2
#I_y = 1
#I_z = 1


Inertia = array( [(I_x   , 0     , 0     ), 
                  (0     , I_y   , 0     ), 
                  (0     , 0     , I_z) ] )

#Inertia = identity(3,float)

#########################################
# SIM
#########################################
t0 = 0
tf = 400 + (1)

dt = 1
ntimes = tf/dt
#print ntimes
x_sim = zeros([0,7])
t = zeros([0])

#print range(t0,ntimes)
for i in range(t0,ntimes):

    if i == 0:
        t_ = range(i*dt, dt*(i+1), dt)
        torque = array([0,0,0.0001])
    else:       
        t_ = range(i*dt-1, dt*(i+1), dt)
        torque = array([0`,0,-0.1*x_sim[-1,2]**2])
        
    #print t_
    #torque_0 = [0.0, 0.0, -0.0001]
    
    
    x_sim_ = integrate.odeint(dynamicsFcn, x_0, t_, args=(torque,Inertia))
    
    t =     concatenate((t,t_),0)
    x_sim = concatenate((x_sim,x_sim_),0)
    #print x_sim
    # for next iteration
    x_0 = x_sim[len(x_sim)-1,:]
    #q_0 = x_sim[-1,3:7]   
    
    #print x_0
    
    #print x_sim
# check the quaternion module
q_mod = sqrt(x_sim[:,3]**2+x_sim[:,4]**2+x_sim[:,5]**2+x_sim[:,6]**2)
#print x_sim

# assign results 
omega_sim = x_sim[:,0:3]
q_sim = x_sim[:,3:7]


# convert all quaternions to euler angles
rpy = zeros((q_sim.shape[0],3))
for i in range(0,q_sim.shape[0]):
    #print q_sim[i,:]
    #rpy[i,:] = tr.euler_from_quaternion(q_sim[i,:])
    rpy[i,:] = quaternion2euler(q_sim[i,:])
    

# assign results 
q_ekf = X[:,0:4]
omega_ekf = X[:,4:7]

#########################################
# INIT FIGURE
#########################################

#turn interactive mode on
ion()

# to get size: fig.get_size_inches()
fig = figure(1) #figsize=(6.5, 9.5))

#clear the figure
fig.clf()
#fig.canvas.manager.window.move(900,0) # in pixels from top-left corner

#ax = fig.add_subplot(111)

#########################################
# PLOTS
#########################################

# plot measurements and compare with estimated fiter measurements
subplot(311)
plot(t, omega_sim*180/pi) # t, omega_meas*180/pi,
ylim(omega_sim.min()*180/pi*1.1-1, omega_sim.max()*180/pi*1.1+1)
legend(['$\omega_x$','$\omega_y$','$\omega_z$'])
ylabel('$\Omega$ [deg/sec]')

subplot(312)
plot(t, q_sim) 
ylim(q_sim.min()*1.1-0.1,q_sim.max()*1.1+0.1)
legend(['$q_1$','$q_2$','$q_3$','$q_4$'])
ylabel('Quaternions')

subplot(313)
plot(t,rpy*180/pi)
#ylim(rpy.min()*180/pi*1.1-1,rpy.max()*180/pi*1.1+1)
legend(['$\phi$ (roll)',r'$\theta$ (pitch)','$\psi$ (yaw)'])
ylabel('Euler Angles [deg]')
#draw()
#show()

#print domega
#print skew(a)

print "elapsed time: {:2}".format((time.time() - start)*1000),"[ms]"
