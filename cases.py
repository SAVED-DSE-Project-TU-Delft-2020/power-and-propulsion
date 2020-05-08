# -*- coding: utf-8 -*-
"""
Created on Thu May  7 12:19:59 2020

@author: Raven
"""

import matplotlib.pyplot as plt
import scipy as np

cd0      = 0.0184
V        = 28
k1       = 0.0245
k2       = 0
g0       = 9.81
WS_range = np.arange(5,200,1)
q        = 457.57

#___CASE 1 CONSTANT ALITITUDE/SPEED CRUISE

def Case1(q, cd0, alpha, WS_range):
    '''
    calculates the wing loading for constant altitude and speed cruise
    inputs are q, cd0, alpha, ws_range
    '''
    T_W = []
    W_S = []
    
    for WS in WS_range:
        TW = (q * cd0) / (alpha * WS)
        T_W.append(TW)
        W_S.append(WS)
    
    plt.plot(W_S, T_W)


#___CASE 2 CONSTANT SPEED CLIMB
        
def Case2(cd0, k1, k2, V, dhdt, q, alpha, beta, WS_range):
    '''
    calculates the wing loading for constant speed climb
    inputs are cd0, k1, k2, cruise V, climb rate dhdt, q, alpha and beta
    '''
    T_W = []
    W_S = []
    
    for WS in WS_range:
        TW =  beta / alpha * (k1 * beta / q * WS + k2 + cd0 / (beta / q * WS) + 1 / V * dhdt)   
        T_W.append(TW)
        W_S.append(WS)
        
    plt.plot(W_S, T_W)
    
    
#___CASE 3 CONSTANT ALTITUDE/SPEED TURN
    
def Case3(cd0, k1, k2, V, g0, Rc, q, alpha, beta, WS_range):
    '''
    calculates the wing loading for constant altitude/ speed turn
    inputs are cd0, k1, k2, cruise V, gravity g0, turn radius Rc, q, alpha and beta
    '''
    T_W = []
    W_S = []
    
    n   = np.sqrt(1+((V**2)/(g0*Rc))**2)
    
    for WS in WS_range:
        TW = beta / alpha * (k1 * n**2 * beta / q * WS + k2 * n + cd0 / (beta / q * WS))
        T_W.append(TW)
        W_S.append(WS)

    plt.plot(W_S, T_W)

    
#___CASE 4 HORIZONTAL ACCELERATION

def Case4(cd0, k1, k2, dVdt, g0, q, alpha, beta, WS_range):
    
    T_W = []
    W_S = []
    
    for WS in WS_range:
        TW = beta / alpha * (k1 * beta / q * WS + k2 + cd0 / beta * q * WS + 1 / g0 * dVdt) 
        T_W.append(TW)
        W_S.append(WS)

    plt.plot(W_S, T_W)

#___CASE 8 SERVICE CEILING
  
def Case8(cd0, k1, k2, V, dhdt, CL, q, alpha, beta, WS_range):

    T_W = []
    W_S = []
    
    for WS in WS_range:
        TW = beta / alpha * (k1 * beta / q * WS + k2 + cd0 / (beta / q * WS) + 1 / V * dhdt)
        T_W.append(TW)
        W_S.append(WS)

    plt.plot(W_S, T_W)
    

#___CASE 0 STALL SPEED
    
def Case0(rho, Vs, CLmax):
    
    WS = 1/2 * rho * Vs**2 * CLmax
    plt.axvline(x=WS)



    
    
    
    
    
    
    
    
    
