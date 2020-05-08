# -*- coding: utf-8 -*-
"""
Created on Fri May  8 15:03:35 2020

@author: Raven
"""

from isacalc import isa
import matplotlib.pyplot as plt
import scipy as np


def Energy1(V, h, cd0, alpha, WS, MTOW, g0, eta1, rang=75*1000,MPL=3):
    
    #dynamic pressure
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    #Calculating power loading
    
    TW   = (q * cd0) / (alpha * WS)
    WP   = 1/(V*TW)
    
    
    #Total efficiency
    
    etatot= eta1
    
    #____DEPARTURE
    
    PreqD = MTOW * g0 / WP
    
    #Available power after efficiencies
    
 
    PbatD  = PreqD * etatot
    
    #Calculating time in cruise and energy consumption
    
    t    = rang/V
    E    = PbatD * t

    #____RETURN
        
    PreqR = (MTOW-MPL) * g0 / WP
    
    #Available power after efficiencies
    
    PbatR  = PreqR * etatot
    
    #Calculating time in cruise and energy consumption
    
    t    = rang/V
    E    += PbatR * t

    return E
    
    
    
    