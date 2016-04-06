# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 10:19:42 2013

@author: miguel
"""

from numpy import *         
from numpy.linalg import * # for inv

#########################################
# Discrete Linear Kalman Filter function
#########################################
# A - State transition matrix
# B - Control matrix
# H - measurement matrx
# X - State vector
# P - State Error Covariance matrix
# Q - Process Error Covariance Matrix
# R - Measuremtn Error Covatiance Matrix
# z - Measurement Vector

def extended_kalman_filter(A, B, H, X, P, Q, R, z, u, dt):
    # EFK discrete and linear
    # [np.newaxis]) is used in case of 1D vectors so they can be transposed
    
    # make sure the vectors are in the right shape
    
    #take care of the 1-D arrays so they are in column shape (this way there are not surprises)
    X = X.reshape(size(X),1)
    z = z.reshape(size(z),1)
    u = u.reshape(size(u),1)
    #H = H[np.newaxis]
    
    ######################################
    ## 1. Prediction Step
    ######################################  
    ## State Prediction (Predict where we're gonna be in the next step)
    # X[n x 1] = A*X + B*u 
    X_   = add(dot(A*dt+eye(size(X)), X), dot(B,u))
    
    ## Covariance Prediction (Predict how much error) 
    # P[n x n] = A*P*A.T + Q  
    P_   = dot(dot(A, P), A.T) + Q

    ######################################
    ## 2. Observation Step
    ######################################   
    
    # Innovation or measurement residual (Compare reality against prediction)
    # y[nz x 1] = z - H*X_
    y = add(z,-dot(H, X_))   
    
    # Innovation Covariance
    # S[nz x nz] = H*P_*H.T + R
    S = dot(dot(H, P_), H.T) + R   
            
    ######################################
    ## 3. Update Step (with Kalman Gain)
    ######################################   
    
    # Kalman Gain (Moderate the prediction)
    # K [n x nz] = P_*H.T*S^(-1)
    # one dimensional case, 1/(...)
    # K = dot(dot(P_[k+1,:,:], H), 1/(dot(dot(H, P_[k+1,:,:]), H.T) + R))
    # multiple dimensional case, inv(...)
    K = dot(dot(P_, H.T), inv(S))
    
    # State Update (New estimate of where we are)
    # X = X_ + K*y
    X = add(X_,dot(K, y)) #[np.newaxis]
    
    # Covariance Update (New estimate of error)
    # one dimensional case, 1/(...)
    # P[k+1,:,:] = dot((identity(3) - K.reshape(size(K),1)*H ), P_[k+1,:,:])
    P = dot((identity(shape(X_)[0]) - dot(K,H) ), P_)
    
    # return the state and covariance estimates
    # X can be reshaped for regular usage in row form
    return X.T.ravel(), P, K