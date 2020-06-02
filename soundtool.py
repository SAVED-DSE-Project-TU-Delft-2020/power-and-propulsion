# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:52:12 2020

@author: Raven
"""

import scipy as np
from scipy import special
from isacalc import isa
import matplotlib.pyplot as plt

#ORDERED ROTATIONAL NOISE

def p_m(m, S, R, P_h, T, B, M_t, theta):
    
    A = R**2*np.pi
    T = 0.22480894244319*T
    J_mb = np.special.jn(m*B, 0.8*M_t*m*B*np.sin(theta))
#    p = 169.3*m*B*R*M_t/(S*A)*(0.76*P_h/M_t**2-np.absolute(T*np.cos(theta)))*J_mb
    
#    p2 = 169.3*m*B*R*M_t/(S*A)*(0.76*P_h/((0.8*M_t)**2)-np.absolute(T*np.cos(theta)))*J_mb
  
    p = 169.3*m*B*R*M_t/(S*A)*(0.76*P_h/M_t**2)*J_mb
    
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

def SPL_to_PWL_dist(SPL, S, S_obs=0):
    
    S = S*0.3048
    A = 4*np.pi*S**2
    if S_obs==0:
        PWL = SPL+10*np.log10(A)
        return PWL
    else:
        S_obs = S_obs*0.3048
        A_obs = 4*np.pi*S_obs**2
        ratio = A/A_obs
        SPL_new = SPL+10*np.log10(ratio)
        return SPL_new


#UNIT CONVERSIONS

def A_ft(d_inch):
    
    d_ft = 0.0833333333* d_inch
    A_ft = d_ft**2*np.pi/4
    
    return A_ft

def V_r(RPM, d_inch, r):
    
    rads = 0.104719755*RPM
    r_inch = d_inch/2
    r_ft = 0.0833333333*r_inch
    
    V_r = rads*r_ft*r
    
    return V_r

#POSITION VECTORS
spanloc_out = 0.7/2
spanloc_in = 0.35/2
b = 3 

engine_right_out = [0, spanloc_out*3, 0]
engine_right_in = [0, spanloc_in*3, 0]
engine_left_out = [0, -spanloc_out*3, 0]
engine_left_in = [0, -spanloc_in*3, 0]
xyzsource = np.array((engine_right_out, engine_right_in, engine_left_out, engine_left_in))
xyzobserver = np.array([-1,0,0])

def position(xyzsource, xyzobserver):
    S = np.array([])
    thetas = np.array([])
    for xyz in xyzsource:
        distance = np.linalg.norm(np.array(xyz)-np.array(xyzobserver))
        opposite = np.absolute(xyz[0]-xyzobserver[0])
        adjacent = np.absolute(xyz[1]-xyzobserver[1])
        angle = np.arctan2(opposite, adjacent)
        S = np.append(S,distance/0.3048)
        thetas = np.append(thetas, angle)
    
    return S, thetas

#NOISE CALCULATION

def noise(xyzsource, xyzobserver, d_inch, P_h, T, B, RPM, h):
    
    S, theta = position(xyzsource, xyzobserver)
    R_ft = 0.0833333333*d_inch/2
    m = 1
    gamma, R, Temp = 1.4, 287, isa(h)[1]
    M_t = V_r(RPM, d_inch, 1)*0.3048/np.sqrt(gamma*R*Temp)
    V_07 = V_r(RPM, d_inch, 0.7) 
    
    #Assumption blade area
    A_b = 2/B*1/9
    
    p_rn = []
    p_v = []
    
    for S, theta in zip(S, theta):
        #ORDERED ROTATIONAL NOISE PRESSURE AT OBSERVER
        rotationp = p_m(m, S, R_ft, P_h, T, B, M_t, theta)
        p_rn.append(rotationp)
        
        #VORTEX NOISE 
        vortexSPL = SPL_vortex(A_b, V_07)
        vortexSPL_obs = SPL_to_PWL_dist(vortexSPL, 300, S_obs=S)
        vortexp = SPL_to_p(vortexSPL_obs)
        p_v.append(vortexp)

        
    p_sum = sum(p_rn+p_v)
    SPL = p_to_SPL(p_sum)

    return SPL



#plst=list()
#thetalst=list()
#for theta in np.arange(0,np.pi*1.05,0.05*np.pi):
#    p=p_m(1, 3, 0.0833333333*15.5/2, 1.36, 55.77, 3, 0.75, theta)
#    thetalst.append(theta)
#    plst.append(p)
#    
#plt.plot(thetalst,plst)
#plt.show()
#plt.ylabel('p')
#plt.xlabel('theta')
#    
    


#plst=list()
#thetalst=list()
#for B in np.arange(0,4.01,1):
#    p,J_mb=p_m(1, 4.757217847769028, 0.645833333075, 1.36, 55.7, B, 0.7053836047832067, 0.7610127542247298)
#    thetalst.append(B)
#    plst.append(B*J_mb)
#    
#plt.plot(thetalst,plst)
#plt.show()
#plt.ylabel('b*bessel')
#plt.xlabel('b')



