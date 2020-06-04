# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:52:12 2020

@author: Raven
"""

import scipy as np
from scipy import special
from isacalc import isa
import matplotlib.pyplot as plt


#POSITIONS
spanloc_out = 0.7/2
spanloc_in = 0.35/2
b = 3 

engine_right_out = [0.4613491583980512, spanloc_out*3, 0]
engine_right_in = [0.2306745791990256, spanloc_in*3, 0]
engine_left_out = [0.4613491583980512, -spanloc_out*3, 0]
engine_left_in = [0.2306745791990256, -spanloc_in*3, 0]
xyzsource = np.array((engine_right_out, engine_right_in, engine_left_out, engine_left_in))
FW_skeleton = [[0,0.6590702262829303,0.9027022644938187,0.6960915377453957,0.9027022644938187,0.6590702262829303,0],[0,-1.5,-1.5,0,1.5,1.5,0]]
CG_loc = [-1,0]


def p_m(m, S, R, P_h, T, B, M_t, theta):
    '''
    Calculates the rotational noise at a certain position (thrust/ power noise)
    '''
    T = 0.22480894244319*T
    theta = theta+0.5*np.pi
    A = R**2*np.pi
    J_mb = np.special.jn(m*B, 0.8*M_t*m*B*np.sin(theta))
    p = 169.3*m*B*R*M_t/(S*A)*(0.76*P_h/M_t**2-T*np.cos(theta))*J_mb
    
    return np.absolute(p)

    
def SPL_vortex(A_b, V_07):
    ''' 
    Calculates the vortex noise at 300 ft
    '''
    k = 3.8*10**(-27)
    SPL = 10*np.log10(k*A_b*V_07**6/(10**-16))
    
    return SPL

#NOISE CONVERSIONS
    
def p_to_SPL(p):
    '''
    Converts a pressure level to SPL
    '''
    p_ref = 0.0002 #dynes/cm^2
    SPL = 20*np.log10(p/p_ref)
    
    return SPL

def SPL_to_p(SPL):
    '''
    Converts SPL to a pressure level
    '''
    p_ref = 0.0002 #dynes/cm^2
    p = 10**(SPL/20)*p_ref
    
    return p

def SPL_to_PWL_dist(SPL, S, S_obs=0):
    '''
    Converts SPL to other distances, or to a power level. S in ft
    '''
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

def PWL_to_SPL(PWL, S):
    '''
    Calculates noise pressure based on noise power source. S in m
    '''
    A = 4*np.pi*S**2
    SPL = PWL-10*np.log10(A)
    
    return SPL
    
#UNIT CONVERSIONS

def A_ft(d_inch):
    '''
    Calculates circle area in ft based on diameter in inch
    '''
    d_ft = 0.0833333333* d_inch
    A_ft = d_ft**2*np.pi/4
    
    return A_ft

def V_r(RPM, d_inch, r):
    '''
    Calculates velocity based on radial position and rpm
    '''
    rads = 0.104719755*RPM
    r_inch = d_inch/2
    r_ft = 0.0833333333*r_inch
    
    V_r = rads*r_ft*r
    
    return V_r


def position(xyzsource, xyzobserver):
    '''
    Calculates the norm and angle between a set of points
    '''
    S = np.array([])
    thetas = np.array([])
    for xyz in xyzsource:
        distance = np.linalg.norm(np.array(xyz)-np.array(xyzobserver))
        opposite = xyz[0]-xyzobserver[0]
        adjacent = xyz[1]-xyzobserver[1]
        angle = np.arctan2(opposite, adjacent)
        S = np.append(S,distance/0.3048)
        thetas = np.append(thetas, angle)
    
    return S, thetas

#NOISE CALCULATION

def noise_heat(xyzsource, d_inch, P_h, T, B, RPM, h):
    '''
    Calculates the SPL on a set of points.
    Also produces heat map of the surroundings
    '''
    R_ft = 0.0833333333*d_inch/2
    m = 1
    gamma, R, Temp = 1.4, 287, isa(h)[1]
    M_t = V_r(RPM, d_inch, 1)*0.3048/np.sqrt(gamma*R*Temp)
    V_07 = V_r(RPM, d_inch, 0.7) 
    
    #Assumption blade area
    A_b = 2/B*1/12
    
    SPLlst = list()
    xrange = np.arange(-0.6,1+0.01,0.01)
    yrange = np.arange(-1.6,1.6+0.01,0.01)
    for x in xrange:
        for y in yrange: 
            xyzobserver=[x,y,0]
            S, theta = position(xyzsource, xyzobserver)
            
            p_rn = list()
            p_v = list()
            for S, theta in zip(S, theta):
                #ORDERED ROTATIONAL NOISE PRESSURE AT OBSERVER
                rotationp = p_m(m, S, R_ft, P_h, T, B, M_t, theta)
                p_rn.append(rotationp)
                #VORTEX NOISE AT OBSERVER
                vortexSPL = SPL_vortex(A_b, V_07)
                vortexSPL_obs = SPL_to_PWL_dist(vortexSPL, 300, S_obs=S)
                vortexp = SPL_to_p(vortexSPL_obs)
                p_v.append(vortexp)


            p_sum = sum(p_rn)+sum(p_v)
            SPL = p_to_SPL(p_sum)
            if SPL>120:
                SPLlst.append(np.nan)
            else:
                SPLlst.append(SPL)

    SPLlst = np.array(SPLlst).reshape(len(xrange),len(yrange))
    window = [yrange[0],yrange[-1],xrange[-1],xrange[0]]
    plt.imshow(SPLlst, cmap='Reds',extent=window)
    plt.colorbar()
    plt.plot(FW_skeleton[1],FW_skeleton[0],"black")
    return SPLlst


def req_noise(xyzsource, d_inch, P_h, T, B, RPM, h, r, CG, PWL_req):
    
    SPL_req = PWL_to_SPL(PWL_req,r)
    print(SPL_req)
    SPLlst = list()
    
    R_ft = 0.0833333333*d_inch/2
    m = 1
    gamma, R, Temp = 1.4, 287, isa(h)[1]
    M_t = V_r(RPM, d_inch, 1)*0.3048/np.sqrt(gamma*R*Temp)
    V_07 = V_r(RPM, d_inch, 0.7) 
    
    #Assumption blade area
    A_b = 2/B*1/12
    
    for theta in np.arange(0,360+1,1):
        
        x = r*np.sin(theta)
        y = r*np.cos(theta)
        
        xyzobserver=[x,y,0]
        S, theta = position(xyzsource, xyzobserver)
    
        p_rn = list()
        p_v = list()
        
        for S, theta in zip(S, theta):
            #ORDERED ROTATIONAL NOISE PRESSURE AT OBSERVER
            rotationp = p_m(m, S, R_ft, P_h, T, B, M_t, theta)
            p_rn.append(rotationp)
            #VORTEX NOISE AT OBSERVER
            vortexSPL = SPL_vortex(A_b, V_07)
            vortexSPL_obs = SPL_to_PWL_dist(vortexSPL, 300, S_obs=S)
            vortexp = SPL_to_p(vortexSPL_obs)
            p_v.append(vortexp)
        
        p_sum = sum(p_rn+p_v)
        
        SPL = p_to_SPL(p_sum)
        SPLlst.append(SPL)
        
    return max(SPLlst)
    
    

#BUTTERFLY PLOT
#plst=list()
#thetalst=list()
#for theta in np.arange(0,2*np.pi*1.026,0.05):
#    p=p_m(1,3,0.645833333075,1.36,12.521858094085683,3, 0.75, theta)
#    if p>0:
#        thetalst.append(theta)
#        plst.append(p)
#    else:
#        thetalst.append(theta)
#        plst.append(-p)
#
#plt.plot(thetalst,plst)
#plt.polar(thetalst,plst)
#plt.show()
#plt.ylabel('p')
#plt.xlabel('theta')