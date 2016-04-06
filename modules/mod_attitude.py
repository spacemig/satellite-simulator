# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 15:46:56 2013

@author: miguelnunes

notes on the attitude code
- for books and papers the quaternion representation with scalar last is common
  but for programming it's more common to find scalar first representation
- to avoid confusion I use: 
    q =     [q1,q2,q3, q4] for scalar last (q4) 
  and 
    q = [q0, q1,q2,q4] for scalar first (q0) 
    
  notice that q1,q2,q3 are the vector components for both cases
  
  references:
      [1] - Analytical Mechanics of Space Systems, Shaub & Junkins, 2003

"""

from numpy import *
from numpy.linalg import *
from scipy import integrate
#import matplotlib.pyplot as plt
#from matplotlib import plot, ion, show
from matplotlib.pylab import *
#import transformations as tr
#from transformations import *

#########################################
# FUNCTIONS
#########################################
#def skew(x):
#    # quaternion representation, scalar last, q=(q1,q2,q3,q4_scalar)
#    return array( [ ( 0   , -x[2],  x[1]), 
#                    ( x[2],  0,    -x[0]),
#                    (-x[1],  x[0],  0   ) ] )

def skew(x):
    # quaternion representation, scalar first, q=(q0,q1,q2,q3)
    # MUST CHECK!!! now it's probvably the same as scalar last!
    return array( [ ( 0   , -x[2],  x[1]), 
                    ( x[2],  0,    -x[0]),
                    (-x[1],  x[0],  0   ) ] )
#def skew4(x):
#    # quaternion representation, scalar last, q=(q1,q2,q3,q4_scalar)
#    # eqn 2.1.7, pg 8
#    return array( [ (  0   ,  x[2], -x[1], x[0]), 
#                    ( -x[2],  0   ,  x[0], x[1]),
#                    (  x[1], -x[0],  0   , x[2]),
#                    ( -x[0], -x[1], -x[2], 0) ] )

def skew4(x):
    # quaternion representation, scalar first, q=(q0,q1,q2,q3)
    # reference [1.100]
    w1 = x[0]
    w2 = x[1]
    w3 = x[2]
    return array( [  [  0,  -w1, -w2, -w3], 
                     [  w1,  0 ,  w3, -w2],
                     [  w2, -w3,  0 ,  w1],
                     [  w3,  w2, -w1,  0 ] ] )

#def euler2quaternion_test(euler):
#    # from http://www.gamedev.net/topic/597324-quaternion-to-euler-angles-and-back-why-is-the-rotation-changing/
#    eX,eY,eZ = euler[0],euler[1],euler[2]
#   
#    c1 = cos(eX/2);
#    s1 = sin(eX/2);
#    
#    c2 = cos(eY/2);
#    s2 = sin(eY/2);
#    
#    c3 = cos(eZ/2);
#    s3 = sin(eZ/2);
#    
#    qx = s1*c2*c3 + c1*s2*s3;
#    qy = c1*s2*c3 - s1*c2*s3;
#    qz = c1*c2*s3 + s1*s2*c3; 
#    
#    qw = c1*c2*c3 - s1*s2*s3;
#    
#    return array([qx,qy,qz,qw])
    
def quaternion_from_euler(euler):
    """
    Convert euler angles to a quaternion with scalar first representation
    """
    # aerospace standard Euler angles sequence Z,Y,X = yaw, pitch, roll
    # 1, psi   - z (yaw,heading)
    # 2, theta - y (pitch)
    # 3, phi   - x (roll)
    
    # check if pitch is not in [-90, 90] deg domain
    if euler[1] >= pi/2:
        print ">>> WARNING! Pitch is more than 90 deg. Results may not be accurate"
        
    if euler[1] <= -pi/2:
        print ">>> WARNING! Pitch is less than -90 deg. Results may not be accurate"
                
        
    #angles = array([r, p, y])
    c = cos( euler/2. )
    s = sin( euler/2. )
    
    # formulat from  Space Vehicle Dynamics and Control, Wie, pg 338
    # q1,q2,q3,q4/scalar
    #q = [s(:,1).*c(:,2).*c(:,3) - c(:,1).*s(:,2).*s(:,3), ...
    # c(:,1).*s(:,2).*c(:,3) + s(:,1).*c(:,2).*s(:,3), ...
    # c(:,1).*c(:,2).*s(:,3) - s(:,1).*s(:,2).*c(:,3), ...
    # c(:,1).*c(:,2).*c(:,3) + s(:,1).*s(:,2).*s(:,3)];
    
    # eqn A.2.15
#    q1 = s[0]*c[1]*c[2] - c[0]*s[1]*s[2]
#    q2 = c[0]*s[1]*c[2] + s[0]*c[1]*s[2]
#    q3 = c[0]*c[1]*s[2] - s[0]*s[1]*c[2]
#    q4 = c[0]*c[1]*c[2] + s[0]*s[1]*s[2]

    # from book: Quaternions and Rotation Sequences pg 167     
    # scalar first representation   
    q0 = c[0]*c[1]*c[2] + s[0]*s[1]*s[2]
    q1 = c[0]*c[1]*s[2] - s[0]*s[1]*c[2]
    q2 = c[0]*s[1]*c[2] + s[0]*c[1]*s[2]
    q3 = s[0]*c[1]*c[2] - c[0]*s[1]*s[2]

    #scalar first
    return array([ q0, q1, q2, q3])


# before quaternion2euler_aero
def euler_from_quaternion(q):
    # from book: Quaternions and Rotation Sequences pg 168

    dcm = dcm_from_quaternion(q)
    
    psi   = arctan2(dcm[0,1],dcm[0,0]) #yaw
    theta = arcsin(-dcm[0,2])         #pitch 
    phi   = arctan2(dcm[1,2],dcm[2,2]) #roll 
    
    return array([psi,theta,phi])
    

#def quaternion2euler_test(q):
#    # from http://www.gamedev.net/topic/597324-quaternion-to-euler-angles-and-back-why-is-the-rotation-changing/
#        
#    qx,qy,qz,qw = q[0],q[1],q[2],q[3]
#    
#    eX = arctan2(-2*(qy*qz-qw*qx), qw*qw-qx*qx-qy*qy+qz*qz)
#    eY = arcsin(2*(qx*qz + qw*qy))
#    eZ = arctan2(-2*(qx*qy-qw*qz), qw*qw+qx*qx-qy*qy-qz*qz)
#    
#    return array([eX,eY,eZ])
    
def quatNorm(q):
    '''
    normalize quaternion
    '''
    return q/sqrt(dot(q,q))
    
# ------------------------------------------------------------------------------
# DCM operations

def dcm_from_euler(euler):
    # from book: Quaternions and Rotation Sequences pg 167
    
    psi, theta, phi = euler
    
    cpsi   = cos(psi)
    spsi   = sin(psi)
    
    ctheta = cos(theta)
    stheta = sin(theta)
    
    cphi   = cos(phi)    
    sphi   = sin(phi)
    
    return array([
    [cpsi*ctheta                  , spsi*ctheta                  , -stheta     ],
    [cpsi*stheta*sphi - spsi*cphi , spsi*stheta*sphi + cpsi*cphi ,  ctheta*sphi],
    [cpsi*stheta*cphi + spsi*sphi , spsi*stheta*cphi - cpsi*sphi ,  ctheta*cphi]
    ])
    

#def dcm_from_quaternion(q):
#    #eqn A.2.13
#    # from Wertz pg 414
#    # q = [q1,q2,q3,q4_scalar)]
#    q1,q2,q3,q4 = q[0],q[1],q[2],q[3]
#    
#    return array([
#                [q1**2 - q2**2 -q3**2 + q4**2,     2*(q1*q2-q3*q4),                 2*(q1*q3+q2*q4)],
#                [2*(q1*q2+q3*q4),                 -q1**2 + q2**2 - q3**2 + q4**2,   2*(q2*q3-q1*q4)],
#                [2*(q1*q3-q2*q4),                  2*(q2*q3+q1*q4),                -q1**2 - q2**2 + q3**2 + q4**2]
#                ])

# before: quaternion2dcm_aero    
def dcm_from_quaternion(q):
    # from book: Quaternions and Rotation Sequences pg 168
    q0,q1,q2,q3 = q #[0],q[1],q[2],q[3]
    
    return array([
                [2*q0**2-1+2*q1**2, 2*(q1*q2+q0*q3),     2*(q1*q3-q0*q2)],
                [2*(q1*q2-q0*q3),   2*q0**2-1+2*q2**2,   2*(q2*q3+q0*q1)],
                [2*(q1*q3+q0*q2),   2*(q2*q3-q0*q1),     2*q0**2-1+2*q3**2]
                ])

###############################################################################
# Rotations

# difference between passive and active transformation
# passive transformation - rotates the frame, but the point/vector remains fixed
# active transformation  - rotates the point/vector but the frame remains fixed
# http://en.wikipedia.org/wiki/Active_and_passive_transformation
def rotX3d_passive(angle):
    # coordinate frame transformation (passive transformation/rotation) 
    # through the given angle 
    cang = cos(angle)
    sang = sin(angle)
    
    return array([
    [1,   0   , 0    ],
    [0,   cang, sang ],
    [0 , -sang, cang ],
    ])
def rotY3d_passive(angle):
    # coordinate frame transformation (passive transformation/rotation) 
    # through the given angle 
    cang = cos(angle)
    sang = sin(angle)
    
    return array([
    [ cang , 0   , -sang],
    [ 0    , 1   , 0   ],
    [ sang , 0   , cang]
    ])
def rotZ3d_passive(angle):
    # coordinate frame transformation (passive transformation/rotation) 
    # through the given angle 
    cang = cos(angle)
    sang = sin(angle)
    
    return array([
    [ cang, sang, 0],
    [-sang, cang, 0],
    [ 0   , 0   , 1]
    ])


def rotX3d_active(angle):
    # point transfromation (active transformation/rotation)
    # through the given angle
    
    # it's just the inverse or transpose of the passive transformation
    return rotX3d_passive(angle).T

def rotY3d_active(angle):
    # point transfromation (active transformation/rotation)
    # through the given angle
    
    # it's just the inverse or transpose of the passive transformation
    return rotY3d_passive(angle).T
    
def rotZ3d_active(angle):
    # point transfromation (active transformation/rotation)
    # through the given angle
    
    # it's just the inverse or transpose of the passive transformation
    return rotZ3d_passive(angle).T


    
#########################################
# DYNAMICS + KINEMATICS
#########################################
# dynamics
#omega  = array( [1,2,3] )
#torque = array( [0,0,0] )
#domega = dot(inv(I), dot( dot(-skew(omega), Inertia), omega) + torque )

def attitude_dynamics(X, t, Torque, Inertia):
    #global Inertia
    
    #state vector x = (q0,q1,q2,q3, wx,wy,wz)
    q     = X[0:4]
    omega_b_i = X[4:7]

    #satellite kinematics
    #omega_b_o is angular velocity of body frame wrt orbital???? or inertial .. frame 
    # represented in body frame - IMU
    #omega_b_i = omega
    dq = 0.5*dot(skew4(omega_b_i),q)

    #satellite dynamics, pg 6
    # inv(skew_omega * Inertia) * omega + torque
    domega = dot(inv(Inertia), dot( dot(-skew(omega_b_i), Inertia), omega_b_i) + Torque)

    
    dX = concatenate((dq,domega),1)
    #print t,dx
    # return state vector q, omega
    return dX

#########################################
# Linearization
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
def linear_dynamics(omega,inertia):
    
    I_x = inertia[0,0]
    I_y = inertia[1,1]
    I_z = inertia[2,2]
    
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