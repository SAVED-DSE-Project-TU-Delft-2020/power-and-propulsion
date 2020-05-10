# -*- coding: utf-8 -*-
"""
Created on Fri May  8 15:03:35 2020

@author: Raven
"""

from isacalc import isa
import matplotlib.pyplot as plt
import scipy as np


def Energy1(V, h, cd0, alpha, WS, MTOW, etaBat, etaElec, DoD, g0=9.81, rang=75*1000,MPL=3):
    '''
    Calculates the total energy consumption in cruise.
    Corresponds to cases 1. Acceleration is neglected
    Returns total energy consumption for return trip
    '''
    #dynamic pressure
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    #Calculating power loading
    
    TW   = (q * cd0) / (alpha * WS)
    WP   = 1/(V*TW)
       
    #Total efficiency and time per run
    
    etatot= etaBat * etaElec
    t    = rang/V
        
    #____DEPARTURE
    
    PreqD = MTOW * g0 / WP
    
    #Available power after efficiencies
    
    PbatD  = PreqD * etatot
    
    #Calculating time in cruise and energy consumption
    

    E    = PbatD * t / DoD

    #____RETURN
        
    PreqR = (MTOW-MPL) * g0 / WP
    
    #Available power after efficiencies
    
    PbatR  = PreqR * etatot
    
    #Calculating time in cruise and energy consumption

    E    += PbatR * t / DoD

    return E



def Energy2(cd0, k1, k2, V, h, alpha, beta, WS, MTOW, etaBat, etaElec, DoD, g0=9.81, MPL=3):
    '''
    Calculates the total energy consumption in VTOL.
    Corresponds to cases 2. Acceleration is negelected
    Returns total energy consumption for return trip
    '''
    #dynamic pressure
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    #Calculating power loading
    TW   =  beta / alpha * (k1 * beta / q * WS + k2 + cd0 / (beta / q * WS) + 1)
    WP   = 1/(V*TW)
    
    #Total efficiency and time per run
    
    etatot= etaBat * etaElec
    t    = h/V
    
    #____DEPARTURE
    
    PreqD = MTOW * g0 / WP
    
    #Available power after efficiencies
    
    PbatD  = PreqD * etatot
    
    #Calculating time in ascend and descend energy consumption
    
    E    = PbatD * t * 2 / DoD

    #____RETURN
        
    PreqR = (MTOW-MPL) * g0 / WP
    
    #Available power after efficiencies
    
    PbatR  = PreqR * etatot
    
    #Calculating time in cruise and energy consumption
    
    E    += PbatR * t * 2 / DoD

    return E
    
    
    
    
    
    
    
    
    
    
    
    