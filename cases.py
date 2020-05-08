# -*- coding: utf-8 -*-
"""
Created on Thu May  7 12:19:59 2020

@author: Raven
"""

from isacalc import isa
import matplotlib.pyplot as plt
import scipy as np

cd0      = 0.0184
k1       = 0.0245
k2       = 0
g0       = 9.81
WS_range = np.arange(5,200,1)


#___CASE 1 CONSTANT ALITITUDE/SPEED CRUISE

def Case1(V, h, cd0, alpha, WS_range):
    '''
    calculates the wing loading for constant altitude and speed cruise
    designing for cruise speed
    '''
    T_W = []
    W_S = []
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    for WS in WS_range:
        TW = (q * cd0) / (alpha * WS)
        T_W.append(1/(V*TW))
        W_S.append(WS)
    
    plt.plot(W_S, T_W)


#___CASE 2 CONSTANT SPEED CLIMB
        
def Case2(cd0, k1, k2, V, h, alpha, beta, WS_range):
    '''
    calculates the wing loading for constant speed climb
    designing for climb 
    '''
    T_W = []
    W_S = []
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    for WS in WS_range:
        TW =  beta / alpha * (k1 * beta / q * WS + k2 + cd0 / (beta / q * WS) + 1)   
        T_W.append(1/(V*TW))
        W_S.append(WS)
        
    plt.plot(W_S, T_W)
    
    
#___CASE 3 CONSTANT ALTITUDE/SPEED TURN
    
def Case3(cd0, k1, k2, V, g0, Rc, h, alpha, beta, WS_range):
    '''
    calculates the wing loading for constant altitude/ speed turn
    designing for turn radius
    '''
    T_W = []
    W_S = []
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    n   = np.sqrt(1+((V**2)/(g0*Rc))**2)
    
    for WS in WS_range:
        TW = beta / alpha * (k1 * n**2 * beta / q * WS + k2 * n + cd0 / (beta / q * WS))
        T_W.append(1/(V*TW))
        W_S.append(WS)

    plt.plot(W_S, T_W)

    
#___CASE 4 HORIZONTAL ACCELERATION

def Case4(cd0, k1, k2, dVdt, g0, V, h, alpha, beta, WS_range):
    '''
    calculates the wing loading for horizontal acceleration
    designing for dVdt (acceleration)
    '''
    T_W = []
    W_S = []
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    for WS in WS_range:
        TW = beta / alpha * (k1 * beta / q * WS + k2 + cd0 / beta * q * WS + 1 / g0 * dVdt) 
        T_W.append(1/(V*TW))
        W_S.append(WS)

    plt.plot(W_S, T_W)

#___CASE 8 SERVICE CEILING
  
def Case8(cd0, k1, k2, V, dhdt, CL, h, alpha, beta, WS_range):
    '''
    calculates the wing loading for service ceiling
    designing for maximum altitude
    '''
    T_W = []
    W_S = []
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    for WS in WS_range:
        TW = beta / alpha * (k1 * beta / q * WS + k2 + cd0 / (beta / q * WS) + 1 / V * dhdt)
        T_W.append(1/(V*TW))
        W_S.append(WS)

    plt.plot(W_S, T_W)
    

#___CASE 0 STALL SPEED
    
def Case0(rho, Vs, CLmax):
    '''
    calculates the wing loading for stall speed
    designing for minimum airspeed
    '''    
    WS = 1/2 * rho * Vs**2 * CLmax
    plt.axvline(x=WS)

    
#___CASE 00 ACCELERATION CLIMB
        
def Case00(cd0, k1, k2, dVdt, V, h, g0, alpha, beta, WS_range):
    '''
    calculates the wing loading for accelerating climb
    designing for climb 
    '''
    T_W = []
    W_S = []
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    n = 1 + dVdt / g0
    
    for WS in WS_range:
        TW =  q / alpha * (cd0 / WS + k1 * (n * beta / q)**2 * WS + k2 * n * beta / q)  + beta / alpha * (1 + 1 / g0 * dVdt) 
        T_W.append(1/(V*TW))
        W_S.append(WS)
        
    plt.plot(W_S, T_W)    
    

    
    
    
