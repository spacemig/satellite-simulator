# test EKF
# this script is to compare with the C++ code for COSMOS simulation

from mod_estimation import *
from mod_attitude import *

import time
start = time.time()

print "\n#######################################\n"

#########################################
# Initial conditions
#########################################

dt = 1
time_array = range(0,1000+1,dt)
#time_array = range(0,1000,1)

# define initial state vector
# convert angles to radians
yaw   = radians(0.) # yaw
pitch = radians(0.) # pitch
roll  = radians(0.) # roll

ypr = array([yaw,pitch,roll])

q_0 = quaternion_from_euler(ypr)
#print q_0

# intial euler angles
omega_x0 = radians(0.)
omega_y0 = radians(0.)
omega_z0 = radians(0.1)

omega_0 = [omega_x0, omega_y0, omega_z0]

#intial state vector: q0,q1,q2,q3, omega_x0,omega_y0,omega_z0,
#x_0 = [0.0, 0.0, 0.2, 1.0, 0.0, 0.0, 0.0] 
x_0 = q_0.tolist() + omega_0

# initial torque
torque_0 = array([0.0, 0.0, 0.0]) 

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

#########################################
# SIM
#########################################

x_sim = integrate.odeint(attitude_dynamics, x_0, time_array, args=(torque_0,Inertia))

# check the quaternion module
q_mod = sqrt(x_sim[:,3]**2+x_sim[:,4]**2+x_sim[:,5]**2+x_sim[:,6]**2)

# assign results 
q_sim = x_sim[:,0:4]
omega_sim = x_sim[:,4:7]


#########################################
# EKF
#########################################
# n iterations for kalman filter
n = size(time_array)+1

# number of states to estimate
# state vector = [q0 q1 q2 q3, w1 w2 w3], scalar first
nstates = 7

# initialize the state vector to contain all values
Xest_can = zeros([n,nstates])
#Xest_ = zeros([n,nstates])

y = zeros([n,nstates])

#initial covariance estimate
Pest_can  = zeros(([n,nstates,nstates]))
#P_ = zeros(([n,nstates,nstates]))
P_0 = identity(nstates)*1e-5

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
   
               
# measurement matrix for perfect star tracker only
# that measures quaternions
H_ST = concatenate( (identity(4), zeros([4,3])), axis=1)

# measurement matrix for rate gyro
# assume it measures directly
H_RG = concatenate( (zeros([3,4]), identity(3)), axis=1)

H = concatenate( (H_ST,H_RG) )

B = zeros([7,4])

# initialization for EKF
Xest_ = array(x_0)


for k in time_array:
    
    # linearize dynamics
    A = attitude_dynamics_linear(omega_0,Inertia)

    # control input
    u = zeros([4,1])

    # 4. innovation error
    # y is the sensor measurement? yes
    # sinclair ST noise/accuracy is about 0.01 deg, from spec sheet
    # euler2quaternion(0.01*pi/180,0,0) ~= randn(4)/15000
    st_noise = randn(4)/15000.
        
    # standard deviation measured on the gyro: 0.00411 rad/sec ~= randn()/100
    rg_noise = randn(3)/100.
    
    #y[k+1,:]
    #k = 0
    
    z = concatenate( (q_sim[k,:] + st_noise, omega_sim[k,:] + rg_noise ) )
    #omega_eric = array([6.44e-8, 1.59e-8, 0.009999])
    #z = array([0.9984, 2.189e-7, 4.831e-8, 0.057, 6.44e-8, 1.59e-8, 0.009999])    
    
    #Xest_ = array(x_0)
    Xest, Pest, K = extended_kalman_filter(A, B, H, Xest_, P_0, Q, R, z, u, dt)
    
    #Xest_can = 
    Xest_can[k+1,:]   = Xest
    Pest_can[k+1,:,:] = Pest

#########################################
# INIT FIGURE
#########################################

#turn interactive mode on
ion()

# to get size: fig.get_size_inches()
fig = figure(1) #figsize=(6.5, 9.5))

#clear the figure
fig.clf()

#########################################
# PLOTS
#########################################

# plot measurements and compare with estimated fiter measurements
subplot(211)
plot(time_array, omega_sim, time_array, ) # t, degrees(omega_ekf
#ylim(omega_ekf.min()*180/pi*1.1-1, omega_ekf.max()*180/pi*1.1+1)
legend(['$\omega_x$','$\omega_y$','$\omega_z$'])
ylabel('$\Omega$ [deg/sec]')

subplot(212)
plot(time_array,q_sim, ) #t,q_mod
#ylim(q_sim.min()*1.1-0.1,q_sim.max()*1.1+0.1)
legend(['$q_0$','$q_1$','$q_2$','$q_3$'])
ylabel('Quaternions')