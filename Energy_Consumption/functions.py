# -*- coding: utf-8 -*-
"""
Created on Fri May  8 11:13:48 2020

@author: Raven
"""
import math
import scipy as np
import matplotlib.pyplot as plt

#ISA

def isa(h):
    #Constants
    
    T0=288.15
    R=287.
    g0=9.80665
    p0=101325.
#    rho0=1.225

    if 0<=h<=11000:

        T=T0-0.0065*h
        p=p0*(T/T0)**(-g0/(-0.0065*R))
        rho=p/(R*T)
          
    elif 11000<h<=20000:

        T=T0-0.0065*11000
        p=22626*math.exp(-g0/(R*T)*(h-11000))
        rho=p/(R*T)

    elif 20000<h<=32000:
        
        T=T0-0.0065*11000+0.001*(h-20000)
        p=5472*(T/216.65)**(-g0/(0.001*R))
        rho=p/(R*T)

    elif 32000<h<=47000:
        
        T=T0-0.0065*11000+0.001*12000+0.0028*(h-32000)
        p=867*(T/228.65)**(-g0/(0.0028*R))
        rho=p/(R*T)
        
    elif 47000<h<=51000:
            
        T=T0-0.0065*11000+0.001*12000+0.0028*15000
        p=111*math.exp(-g0/(R*T)*(h-47000))
        rho=p/(R*T)
    
    elif 51000<h<=71000:                          
                            
        T=T0-0.0065*11000+0.001*12000+0.0028*15000-0.0028*(h-51000)
        p=67*(T/270.65)**(-g0/(-0.0028*R))
        rho=p/(R*T)

    elif 71000<h<=84852:
                      
                            
        T=T0-0.0065*11000+0.001*12000+0.0028*15000-0.0028*20000-0.002*(h-71000)
        p=67*(T/214.65)**(-g0/(-0.0028*R))
        rho=p/(R*T)

    return p, T, rho


#FLIGHT PHASES

def P_to(W, rho, A_prop, V_to):
    T_to = 1.2 * W
    
    P_to = 0.5*T_to*V_to*np.sqrt(1+(2*T_to/(rho*V_to**2*A_prop)))
    return P_to

def P_climb(W, Cl_max, Cd_climb, S, phi, rho):
    phi = phi * np.pi/180
    rho_sl=1.225
    V_climb = 1.2*np.sqrt(2*W/(S*rho_sl*Cl_max))
    V_ver = V_climb*np.sin(phi)
    q_c = q(rho, V_climb)
    
    P_climb = (W * np.sin(phi) + Cd_climb*q_c*S)*V_climb
    return P_climb, V_ver

def P_cruise(q, S, Cd, V_cruise):
    T_cruise = q*S*Cd
    
    P_cruise = V_cruise*T_cruise
    return P_cruise

def P_landing(W, A_prop, V_des=4):
    rho_sl = 1.225
    K = 1
    V_h = np.sqrt(W/(2*rho_sl*A_prop))
    x = -V_des/V_h
    V_i = (K-1.125*x -1.372*x**2-1.718*x**3-0.655*x**4)*V_h
    
    P_landing = K*W*(V_i-V_des)
    return P_landing

#FLIGHT CONDITIONS
    
def q(rho, V):
    
    q = 0.5*rho*V**2
    return q

def A_prop(d_prop):
    n_prop = 4
    
    A_prop = np.pi()*d_prop**2/4*n_prop
    return A_prop

def V_cruise(W, Cl, S, rho): 
    V_c = np.sqrt(2 * W / (Cl * S * rho)) 
    return V_c

def Cd(Cl, A, e, Cd0):
    
    Cd = Cd0+Cl**2/(np.pi*A*e)
    return Cd

#WEIGHT
    
def W_tot(m_bat, m_eng, m_struc, m_sensors, PL=True):
    if PL:
        m_pl = 3
    else: 
        m_pl = 0
    g = 9.81
    
    W_tot = (m_bat+m_eng+m_struc+m_sensors+m_pl)*g
    return W_tot

#ENERGIES
    
def E_to(W, A_prop, V_to=6, h_trans=20):
    
    h = np.array([0])
    t = np.array([0])
    dt = 0.01
    E = 0

    while h[-1]<h_trans:
        _,_,rho = isa(round(h[-1],4))
        P = P_to(W, rho, A_prop, V_to)
        E += P*dt
        h = np.append(h,h[-1]+V_to*dt)
        t = np.append(t,t[-1]+dt)
   
    return E, t[-1], h[-1]

def E_climb(W, Cl_max, Cd_climb, S, phi, h_cruise=500, h_trans=20):
    
    h = np.array([h_trans])
    t = np.array([0])
    P_c = np.array([0])
    dt = 0.01
    E = 0

    
    while h[-1]<h_cruise:
        _,_,rho = isa(round(h[-1],4))
        P, V_ver = P_climb(W, Cl_max, Cd_climb, S, phi, rho)
        E += P*dt
        h = np.append(h,h[-1]+V_ver*dt)
        t = np.append(t,t[-1]+dt)
        P_c = np.append(P_c,P)
        
    return E, t[-1], h[-1]
        
def E_cruise(W, S, Cl_cruise, LD, A, e, Cd0, h_cruise=500,h_trans=20,phi=45, r=75000):
    phi = phi * np.pi/180
    s_climb = (h_cruise-h_trans)/np.sin(phi)
    s_glide = LD*(h_cruise-h_trans)
    s = r-s_climb-s_glide
    
    x = np.array([0])
    t = np.array([0])
    P_c = np.array([0])
    dt = 0.1
    E = 0
    
    while x[-1]<s:
        _,_,rho = isa(round(h_cruise,4))
        V_c = V_cruise(W, Cl_cruise, S, rho)
        Cd_c = Cd(Cl_cruise, A, e, Cd0)
        q_c = q(rho, V_c) 
        P = P_cruise(q_c, S, Cd_c, V_c)
        E += P*dt
        x = np.append(x,x[-1]+V_c*dt)
        t = np.append(t,t[-1]+dt)
        P_c = np.append(P_c,P)
        
    return E, t[-1]


def E_landing(W, A_prop, h_landing=20, V_des=4):
    
    h = np.array([h_landing])
    t = np.array([0])
    dt = 0.01
    E = 0

    while h[-1]>0:
        _,_,rho = isa(round(h[-1],4))
        P = P_landing(W, A_prop,V_des=4)
        E += P*dt
        h = np.append(h,h[-1]-V_des*dt)
        t = np.append(t,t[-1]+dt)
   
    return E, t[-1], h[-1]
        
    
def E_sensors(t, P_sensors=25.25):
    
    E_s = t*P_sensors
    
    return E_s 


def E_trip(W, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, A, e, Cd0):
    
    E_T,t_t,_ = E_to(W, A_prop, V_to=6, h_trans=20)
    E_Cl,t_Cl,_ = E_climb(W, Cl_max, Cd_climb, S, phi, h_cruise=500, h_trans=20)
    E_Cr,t_Cr = E_cruise(W, S, Cl_cruise, LD, A, e, Cd0, h_cruise=500,h_trans=20,phi=45, r=75000)
    E_L,t_L,_ = E_landing(W, A_prop, h_landing=20, V_des=4)
    
    t = t_t+t_Cl+t_Cr+t_L
    E_s= E_sensors(t)
    
    E = E_T+E_Cl+E_Cr+E_L+E_s   

    return E, t






        
