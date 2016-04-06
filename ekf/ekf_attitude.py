# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 10:31:17 2013

@author: miguelnunes

references: 
    
    MS THesis, Aykut Kutlu
    Design Of Kalman Filter Based Attitude Determination Algorithms For A 
    Leo Satellite And For A Satellite Attitude Control Test Setup

assumptions: quaternion scalar is last vector entry

notes: 
    for quaternion representation makes sense to use q0,q1,q2,q3 for computational code because of the indices (starting with 0)
    for mathematical formulation, etc. makes more sense to use the other representation with scalar last
   
TBD
- 3d representation of the atittude
- add disturbance torques
- implement the sensors
- ekf     
"""


from mod_estimation import *
from mod_attitude import *

import time
start = time.time()

print "\n#######################################\n"

#########################################
# Initial conditions
#########################################

t = range(0,1000,1)

# define initial state vector
# convert angles to radians
yaw   = radians(0.) # yaw
pitch = radians(0.) # pitch
roll  = radians(0.) # roll

ypr = array([yaw,pitch,roll])

#q_0 = tr.quaternion_from_euler(roll,pitch,yaw)
q_0 = quaternion_from_euler(ypr)
#q_0 = array([0,0,0,1])
omega_x0 = 1*pi/180
omega_y0 = 0*pi/180
omega_z0 = 0*pi/180

omega_0 = [omega_x0, omega_y0, omega_z0]
#x_0 = [0.0, 0.0, 0.2, 1.0, 0.0, 0.0, 0.0]  #intial state vector: omega_x0,omega_y0,omega_z0,q0,q1,q2,q3  
x_0 = omega_0 + q_0.tolist()

torque_0 = [0.0, 0.0, 0.0] # initial torque

#########################################
#a = array ([1,2,3])
#B = array( [ (1,2,3), (4,5,6), (7,8,9) ] )

Inertia = array( [(2.4388, 0.0573, -0.0396), 
            (0.0573, 2.4155, -0.0260), 
            (-0.0396, -0.0260, 2.8271) ] )

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
#ID = identity(3,float)

#########################################
# SIM
#########################################

x_sim = integrate.odeint(dynamicsFcn, x_0, t, args=(torque_0,Inertia))

# check the quaternion module
q_mod = sqrt(x_sim[:,3]**2+x_sim[:,4]**2+x_sim[:,5]**2+x_sim[:,6]**2)
#print x_sim

# assign results 
omega_sim = x_sim[:,0:3]
q_sim = x_sim[:,3:7]

#########################################
# EKF
#########################################

# moved to mod_attitude -> attitude_dynamics_linear
# linearization for quaternion scalar last formulation (default!)
def linearSystem(omega):
    
    omega_x = omega[0]
    omega_y = omega[1]
    omega_z = omega[2]
    
    # attitude section
    A_q = 0.5*array([ ( 0      ,  omega_z, -omega_y,  omega_x ),
                      (-omega_z,  0      ,  omega_x,  omega_y ),
                      ( omega_y, -omega_x,  0      ,  omega_z ),
                      (-omega_x, -omega_y, -omega_z,  0       )])
    
    # append the zeros matrix for the attitude part
    A_q = concatenate((A_q,zeros([4,3])),axis=1)


    # dynamics section
    
    k_1 = -(I_x - I_z)/I_y
    k_2 = (I_x - I_y)/I_z
    k_3 = (I_y - I_z)/I_x
    
    A_omega = array([ ( 0          , k_3*omega_z , k_3*omega_y), 
                      ( k_1*omega_z, 0           , k_1*omega_x), 
                      ( k_2*omega_y, k_2*omega_x , 0          ) ]) 

    #print A_omega    
    # append the zeros matrix for the attitude part
    A_omega = concatenate( (zeros([3,4]),A_omega), axis=1)
    
    A = concatenate ((A_q,A_omega), axis=0)
    #print "#############################"
    #print A_omega
    return A



#initial estimate
# estimate only quaternions for now!
# using scalar first

#timeRange = range(0,5,1)

# n iterations for kalman filter
n = size(t)

# number of states to estimate
# state vector = [q1 q2 q3 q4, w1 w2 w3], scalar last
nstates = 7

# initialize the state vector to contain all values
X = zeros([n,nstates])
X_ = zeros([n,nstates])

y = zeros([n,nstates])

#timek = zeros([n])


#initial covariance estimate
P  = zeros(([n,nstates,nstates]))
P_ = zeros(([n,nstates,nstates]))
P[0,:,:] = identity(nstates)*1e-5

# initial covariance matrices
# model covariance, usually 1000x smaller than P0
# smaller Q will make the EKF trust more the linearized model than the sensor
Q_q = identity(4)*1e-10
Q_omega = identity(3)*5e-7

Q = concatenate ( (concatenate( (Q_q, zeros([4,3])), axis=1), 
                  concatenate( (zeros([3,4]), Q_omega), axis=1)))
# sensor covariance
# sinclair ST noise/accuracy is about 0.01 deg 0.0001745 rad, from spec sheet
R_ST = identity(4)*(0.0001745)**2

# standard deviation measured on the gyro: 0.00411 rad/sec 
R_RG = identity(3)*(0.00411)**2

# realistic r
R = concatenate ( (concatenate( (R_ST, zeros([4,3])), axis=1), 
                  concatenate( (zeros([3,4]), R_RG), axis=1)))
# 
#R = identity(nstates)*1

print R
#print X,P


# measureent matrix for perfect star tracke only
# that measures quaternions
H_ST = concatenate( (identity(4), zeros([4,3])), axis=1)

# measurement matrix for rate gyro
# assume it measures directly
H_RG = concatenate( (zeros([3,4]), identity(3)), axis=1)

H = concatenate( (H_ST,H_RG) )

# initial w measured 
# in the practice the initial values should come from the sensors with added noise
#w = omega_0 + array([0.1,0,0])
# initial state vector estimate (quaternion)
#X = array([1,0,0,0])
#X[0,:] = array([1,0,0,0])
X[0,0:4] = q_0 + array([0.1,-1,0,0.2])
X[0,4:7] = omega_0 + array([5*pi/180,2*pi/180,-1*pi/180])
#print X[0,:]

dt = 1 #sec
#print linearSystem(w)
# EKF, pg 66
for k in (t):
    
    if (k == max(t)):
        break
    #print k
    #timek[k+1] = k+1
    
    #####
    # State Transition Model    
    # print A,X
    # linear system
    A = linearSystem(X[k,4:7])
    #A = linearSystem(w)
    #print A
    #####
    # Predict
    # predict state propagation + #normalize quaterion propagation
    
    # 1. Stave vector propagation
    Phi = (identity(nstates) + A)*dt
    X_[k+1,:] = quatNorm(dot(Phi,X[k,:])) # + randn(nstates)/400# add noise later
    
    #print X_[k+1,:]

    # 2. covariance matrix propagation
    P_[k+1,:,:] =  dot(dot(Phi,P[k,:,:]),Phi.transpose()) + Q
    
    # 3. kalman gain
    invH = inv(dot(dot(H,P_[k,:,:]),H.T) + R)
    K = dot(dot(P_[k,:,:],H.T),invH)
    
    # 4. innovation error
    # y is the sensor measurement? yes
    # sinclair ST noise/accuracy is about 0.01 deg, from spec sheet
    # euler2quaternion(0.01*pi/180,0,0) ~= randn(4)/15000
    st_noise = randn(4)/15000
    
    # standard deviation measured on the gyro: 0.00411 rad/sec ~= randn()/100
    rg_noise = randn(3)/100
    y[k+1,:] = concatenate( (q_sim[k,:] + st_noise, omega_sim[k,:] + rg_noise ) ) 
    
    #simulate ST senssor failure
    if (k >= 300 and k <= 550):
        y[k,0:4] = st_noise*1000
    
    e = y[k,:] - dot(H,X_[k+1,:])
    
    # 5. state update    
    #if (k == 60 ):
    #    X_[k+1,:] = concatenate( (q_sim[k,:] , omega_sim[k,:] ) ) 
        
    X[k+1,:] = X_[k+1,:] + dot(K,e)
    
    # 6. covariance update
    P[k+1,:,:] = dot(identity(nstates) - dot(K,H),P_[k+1,:,:])
    #print e   
    
    #print X
    
# state vector propagation
#xhat_ = phi* 

#print X

# check difference between estimated state and real!
error = sqrt( (X[:,0]-q_sim[:,0])**2 + (X[:,1]-q_sim[:,1])**2 + + (X[:,2]-q_sim[:,2])**2 + + (X[:,3]-q_sim[:,3])**2  )
print error


# convert all quaternions to euler angles
rpy = zeros((q_sim.shape[0],3))
for i in range(0,q_sim.shape[0]):
    #print q_sim[i,:]
    #rpy[i,:] = tr.euler_from_quaternion(q_sim[i,:])
    rpy[i,:] = quaternion2euler(q_sim[i,:])
    
# convert all quaternions to euler angles
rpy_X = zeros((X.shape[0],3))
for i in range(0,X.shape[0]):
    #print q_sim[i,:]
    #rpy[i,:] = tr.euler_from_quaternion(q_sim[i,:])
    rpy_X[i,:] = quaternion2euler(X[i,:])
    


# assign results 
q_ekf = X[:,0:4]
omega_ekf = X[:,4:7]

q_meas = y[:,0:4]
omega_meas = y[:,4:7]

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
subplot(411)
plot(t, omega_sim*180/pi, t, omega_ekf*180/pi) # t, omega_meas*180/pi,
ylim(omega_ekf.min()*180/pi*1.1-1, omega_ekf.max()*180/pi*1.1+1)
legend(['$\omega_x$','$\omega_y$','$\omega_z$'])
ylabel('$\Omega$ [deg/sec]')

'''
subplot(512)
plot(t,q_sim) #t,q_mod
ylim(q_sim.min()*1.1-0.1,q_sim.max()*1.1+0.1)
legend(['$q_1$','$q_2$','$q_3$','$q_4$'])
ylabel('Quaternions')
'''
subplot(412)
plot(t, q_meas, t, q_ekf) #t, q_meas
ylim(q_sim.min()*1.1-0.1,q_sim.max()*1.1+0.1)
legend(['$q_1$','$q_2$','$q_3$','$q_4$'])
ylabel('Est. Quaternions')

subplot(413)
plot(t,error)
ylim(-1.1,1.1)
#legend(['$q_1$','$q_2$','$q_3$','$q_4$'])
ylabel('Error. Est.')

subplot(414)
plot(t,rpy*180/pi,t,rpy_X*180/pi)
ylim(rpy.min()*180/pi*1.1-1,rpy.max()*180/pi*1.1+1)
legend(['$\phi$ (roll)',r'$\theta$ (pitch)','$\psi$ (yaw)'])
ylabel('Euler Angles [deg]')
#draw()
#show()

#print domega
#print skew(a)

print "elapsed time: {:2}".format((time.time() - start)*1000),"[ms]"
