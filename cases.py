# -*- coding: utf-8 -*-
"""
Created on Thu May  7 12:19:59 2020

@author: Raven
"""

import matplotlib.pyplot as plt
import scipy as np



#___CASE 1 CONSTANT ALITITUDE/SPEED CRUISE

def Case1(q, cd0, alpha, WS_range):
    
    T_W = []
    W_S = []
    
    for WS in WS_range:
        TW = (q * cd0) / (alpha * WS)
        T_W.append(TW)
        W_S.append(WS)
    
    plt.plot(W_S, T_W)


#___CASE 2 CONSTANT SPEED CLIMB
        
def Case2(cd0, k1, k2, V, dhdt, q, alpha, beta, WS_range):
    
    T_W = []
    W_S = []
    
    for WS in WS_range:
        TW =  beta / alpha * (k1 * beta / q * WS + k2 + cd0 / (beta / q * WS) + 1 / V * dhdt)   
        T_W.append(TW)
        W_S.append(WS)
        
    plt.plot(W_S, T_W)
    
    
#___CASE 3 CONSTANT ALTITUDE/SPEED TURN
    
def Case3(cd0, k1, k2, V, g0, Rc, q, alpha, beta, WS_range):
    
    T_W = []
    W_S = []
    
    n   = np.sqrt(1+((V**2)/(g0*Rc))**2)
    
    for WS in WS_range:
        TW = beta / alpha * (k1 * n**2 * beta / q * WS + k2 * n + cd0 / (beta / q * WS))
        T_W.append(TW)
        W_S.append(WS)

    plt.plot(W_S, T_W)

    
#___CASE 4 HORIZONTAL ACCELERATION

def Case4():
    
    return



#___CASE 8 SERVICE CEILING
  
def Case8(cd0, k1, k2, V, dhdt, CL, q, alpha, beta, WS_range):

    T_W = []
    W_S = []
    
    for WS in WS_range:
        TW = beta / alpha * (k1 * beta / q * WS + k2 + cd0 / (beta / q * WS) + 1 / V * dhdt)
        T_W.append(TW)
        W_S.append(WS)

    plt.plot(W_S, T_W)



    
    
    
    
    
    
    
    
    
