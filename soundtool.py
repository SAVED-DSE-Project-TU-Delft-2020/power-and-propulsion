# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:52:12 2020

@author: Raven
"""

import scipy as np
from scipy import special

#ORDERED ROTATIONAL NOISE

def p_m(m, S, R, A, P_h, T, B, M_t, theta):
    
    J_mb = np.special.jn(m*B, 0.8*M_t*m*B*np.sin(theta))
    p = 169.3*m*B*R*M_t/(S*A)*(0.76*P_h/M_t**2-T*np.cos(theta))*J_mb
    
    return p

#VORTEX NOISE
    
def SPL_vortex(A_b, V_07):
    
    k = 3.8*10**(-27)
    SPL = 10*np.log10(k*A_b*V_07**6/(10**-16))
    
    return SPL

#NOISE CONVERSIONS
    
def p_to_SPL(p):
    
    p_ref = 0.0002 #dynes/cm^2
    SPL = 20*np.log10(p/p_ref)
    
    return SPL

def SPL_to_p(SPL):
    
    p_ref = 0.0002 #dynes/cm^2
    p = 10**(SPL/20)*p_ref
    
    return p

def SPL_to_PWL(SPL,S):
    
    S = S*0.3048
    A = 4*np.pi*S**2
    PWL = SPL+10*np.log10(A)
    
    return PWL

#UNIT CONVERSIONS

def A_ft(d_inch):
    
    d_ft = 0.0833333333* d_inch
    A_ft = d_ft**2*np.pi/4
    
    return A_ft

def V_r(RPM, d_inch):
    
    rads = 0.104719755*RPM
    r_inch = d_inch/2
    r_ft = 0.0833333333*r_inch
    
    V_r = rads*r_ft*0.7
    
    return V_r










    
    
    
    
    
    
    