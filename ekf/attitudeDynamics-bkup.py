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


from ekf_support_lib import *
import time
start = time.time()

print "\n#######################################\n"

#########################################
# Initial conditions
#########################################

t = range(0,100,1)

# define initial state vector
r = 10*pi/180 # roll in deg
p = 0*pi/180 # pithc in deg
y = 0*pi/180 # yaw in deg

#q0 = tr.quaternion_from_euler(r,p,y)
q0 = euler2quaternion(r,p,y)
#q0 = array([0,0,0,1])
wx = 10*pi/180
wy = 1*pi/180
wz = -1*pi/180

w0 = [wx, wy, wz]
#x0 = [0.0, 0.0, 0.2, 1.0, 0.0, 0.0, 0.0]  #intial state vector: wx,wy,wz,q0,q1,q2,q3  
x0 = w0 + q0.tolist()

torque = [0.0, 0.0, 0.0] # initial torque

#########################################
#a = array ([1,2,3])
#B = array( [ (1,2,3), (4,5,6), (7,8,9) ] )

Inertia = array( [(2.4388, 0.0573, -0.0396), 
            (0.0573, 2.4155, -0.0260), 
            (-0.0396, -0.0260, 2.8271) ] )

#using principal moments of inertia
Inertia = array( [(2.4388, 0     , 0     ), 
                  (0     , 2.4155, 0     ), 
                  (0     , 0     , 2.8271) ] )

#Inertia = identity(3,float)
ID = identity(3,float)


#########################################
# SIM
#########################################

x = integrate.odeint(dynamicsFcn, x0, t, args=(torque,Inertia))

# check the quaternion module
q_mod = sqrt(x[:,3]**2+x[:,4]**2+x[:,5]**2+x[:,6]**2)
#print x

# assign results 
omega = x[:,0:3]
q = x[:,3:7]


    

#########################################
# EKF
#########################################


# linearization for quaternion scalar first formulation
def linearSystem_(w):
    #A = zeros([4,4])
    
    '''
    A[0,0] = 0
    A[0,1] = w[0]
    A[0,2] = w[1]
    A[0,3] = w[2]
    
    A[1,0] = w[0]
    A[1,1] = 0
    A[1,2] = w[2]
    A[1,3] = -w[1]
    
    A[2,0] = w[1]
    A[2,1] = -w[2]
    A[2,2] = 0
    A[2,3] = w[0]
    
    A[3,0] = w[2] 
    A[3,1] = w[1]
    A[3,2] = -w[0]
    A[3,3] = 0
    '''
    
    A = 0.5*array([ ( 0   , -w[0], -w[1], -w[2] ),
                    ( w[0],  0   ,  w[2], -w[1] ),
                    ( w[1], -w[2],  0   ,  w[0] ),
                    ( w[2],  w[1], -w[0],  0    )])
     
    #A = 0.5*A
    return A


# linearization for quaternion scalar last formulation (default!)
def linearSystem(w):
    
    A = 0.5*array([ ( 0   ,  w[2], -w[1],  w[0] ),
                    (-w[2],  0   ,  w[0],  w[1] ),
                    ( w[1], -w[0],  0   ,  w[2] ),
                    (-w[0], -w[1], -w[2],  0    )])
     
    return A


def quatNorm(q):
    return q/sqrt(dot(q,q))
    
#initial estimate
# estimate only quaternions for now!
# using scalar first

#timeRange = range(0,5,1)

n = size(t)

# initialize the state vector to contain all values
X = zeros([n,4])
X_ = zeros([n,4])
#timek = zeros([n])

#initial covariance estimate
P = zeros(([n,4,4]))
P_ = zeros(([n,4,4]))

Q = identity(4)/1
R = identity(4)/10

# initial quaternion estimate
#X = array([1,0,0,0])
#X[0,:] = array([1,0,0,0])
X[0,:] = q0 + array([0.1,0,0,0.2])
#print X

#print X,P


# measureent matrix for perfect star tracke only
# that measures quaternions
H = identity(4)

#initial w measured?
w = array([wx,wy,wz]) + array([2.1,0,0])
#w = array([[1],[0],[0]])

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
    #print A,X
    # linear system
    A = linearSystem(w)
    
    #####
    # Predict
    # predict state propagation + #normalize quaterion propagation
    
    # 1. Stave vector propagation
    Phi = (identity(4) + A)*dt
    X_[k+1,:] = quatNorm(dot(Phi,X[k,:])) + randn(4)/40# add noise later
    
    # 2. covariance matrix propagation
    P_[k+1,:] =  dot(dot(Phi,P[k,1]),Phi.transpose()) + Q
    
    # 3. kalman gain
    invH = inv(dot(dot(H,P_[k,:]),H.T)+R)
    K = dot(dot(P_[k,:],H.T),invH)
    
    # 4. innovation error
    # y is the sensor measurement? yes
    y = q[k,:] + randn(4)/30
    
    #simulate sensor failure
    if (k >= 10 and k <= 20):
        y = randn(4)/10
    
    e = y - dot(H,X_[k+1,:])
    
    # 5. state update    
    X[k+1,:] = X_[k+1,:] + dot(K,e)
    
    # 6. covariance update
    P[k+1,:] = dot(identity(4) - dot(K,H),P_[k+1,:])
    #print e   
    
    #print X
    
# state vector propagation
#xhat_ = phi* 

#print X

# check difference between estimated state and real!
error = sqrt( (X[:,0]-q[:,0])**2 + (X[:,1]-q[:,1])**2 + + (X[:,2]-q[:,2])**2 + + (X[:,3]-q[:,3])**2  )
print error


# convert all quaternions to euler angles
rpy = zeros((q.shape[0],3))
for i in range(0,q.shape[0]):
    #print q[i,:]
    #rpy[i,:] = tr.euler_from_quaternion(q[i,:])
    rpy[i,:] = quaternion2euler(q[i,:])
    
# convert all quaternions to euler angles
rpy_X = zeros((X.shape[0],3))
for i in range(0,X.shape[0]):
    #print q[i,:]
    #rpy[i,:] = tr.euler_from_quaternion(q[i,:])
    rpy_X[i,:] = quaternion2euler(X[i,:])
    

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
subplot(411)
plot(t,omega*180/pi)
ylim(omega.min()*180/pi*1.1-1, omega.max()*180/pi*1.1+1)
legend(['$\omega_x$','$\omega_y$','$\omega_z$'])
ylabel('$\Omega$ [deg/sec]')

'''
subplot(512)
plot(t,q) #t,q_mod
ylim(q.min()*1.1-0.1,q.max()*1.1+0.1)
legend(['$q_1$','$q_2$','$q_3$','$q_4$'])
ylabel('Quaternions')
'''
subplot(412)
plot(t,X,t,q)
ylim(q.min()*1.1-0.1,q.max()*1.1+0.1)
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
