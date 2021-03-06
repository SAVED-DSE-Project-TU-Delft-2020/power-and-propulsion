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
    '''
    Calculates pressure, temperature and density based on ISA at certain height.
    '''
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
    '''
    Calculates power for takeoff
    '''
    T_to = 1.2 * W
    P_to = 0.5*T_to*V_to*np.sqrt(1+(2*T_to/(rho*V_to**2*A_prop)))
#    print(P_to, T_to**(3/2)/np.sqrt(2*rho*A_prop))  #Verification turkish paper vs italian paper, Design and prototyping high performance multirotor 
    return P_to, T_to

def P_climb(W, Cl_max, Cd_climb, S, phi, rho):
    '''
    Calculates power for climb
    '''
    phi = phi * np.pi/180
    rho_sl = 1.225
    V_climb = 1.2*np.sqrt(2*W/(S*rho_sl*Cl_max))
    V_ver = V_climb*np.sin(phi)
    q_c = q(rho, V_climb)
    T_climb = W * np.sin(phi) + Cd_climb*q_c*S
    P_climb = T_climb*V_climb
    return P_climb, V_ver, V_climb, T_climb

def P_cruise(q, S, Cd, V_cruise):
    '''
    Calculates power for cruise
    '''
    T_cruise = q*S*Cd
    P_cruise = V_cruise*T_cruise
    return P_cruise, T_cruise

def P_landing(W, A_prop, V_des=4):
    '''
    Calculates power for landing
    '''
    rho_sl = 1.225
    K = 1
    V_h = np.sqrt(W/(2*rho_sl*A_prop))
    x = -V_des/V_h
    V_i = (K-1.125*x -1.372*x**2-1.718*x**3-0.655*x**4)*V_h
    P_landing = K*W*(V_i-V_des)
#    print(P_landing, W**(3/2)/(2*rho_sl*A_prop)) #Verification turkish paper vs italian paper Design and prototyping high performance multirotor, Design and prototyping high performance multirotor
    return P_landing, W

#FLIGHT CONDITIONS
    
def q(rho, V):
    '''
    Calculates dynamic pressure
    '''
    q = 0.5*rho*V**2
    return q

def A_pr(d_prop):
    '''
    Calculates propeller area
    '''
    n_prop = 4
    
    A_pr = np.pi()*d_prop**2/4*n_prop
    return A_pr

def V_cruise(W, Cl, S, rho): 
    '''
    Calculates cruise speed based on L=W
    '''
    V_c = np.sqrt(2 * W / (Cl * S * rho)) 
    return V_c

def Cd(Cl, A, e, Cd0):
    '''
    Calculates Cd using drag polar
    '''
    Cd = Cd0+Cl**2/(np.pi*A*e)
    return Cd

#PROPELLER CONDITIONS
    
def V_r(RPM, d_inch, r):
    '''
    Calculates velocity based on radial position and rpm
    '''
    rads = 0.104719755*RPM
    r_inch = d_inch/2
    r_m = 0.0254*r_inch
    V_r = rads*r_m*r
    return V_r

def v_axial_propeller(V_0,T,rho,A_prop):
    '''
    Calculates the axial velocity
    '''
    A_prop = A_prop/4 #single propeller area
    v_a = 0.5*(-V_0 + np.sqrt(V_0**2 + 2*T/(rho*A_prop))) #Equation A.29 PhD Veldhuis
    return v_a 

def n_prop(V_a,V_cruise): 
    '''
    Calculates the efficiency of the propulsion
    '''
    return  V_cruise/(V_cruise+V_a)

def CT(T,rho,A_prop,V_tip):
    '''
    Calculates the thrust coefficient
    '''
    return T/(rho*A_prop*V_tip**2)

def FM(CT, Cd0, sigma=0.054, k=1.15):
    '''
    Calculates the figure of merit
    '''
    nominator = CT**(3/2)/np.sqrt(2)
    denominator = k*CT**(3/2)/np.sqrt(2)+sigma*Cd0/8
    return nominator/denominator

def RPM(T):
    '''
    Calculates RPM through empirical relation
    '''
    return 106*T+2831
#WEIGHT
    
def W_tot(m_bat, m_eng, m_struc, m_sensors, PL=True):
    '''
    Calculates total weight
    '''
    if PL:
        m_pl = 3
    else: 
        m_pl = 0
    g = 9.81
    
    W_tot = (m_bat+m_eng+m_struc+m_sensors+m_pl)*g
    return W_tot

#ENERGIES
    
def E_to(W, A_prop, Cd0, t=np.array([0]), E=np.array([0]), P_arr=np.array([0]), V_to=6, h_trans=20):
    '''
    Calculates energy for takeoff
    '''
    t_begin = t[-1]
    h = np.array([0.])
    dt = 0.01
    
    while h[-1]<h_trans:
        _,_,rho = isa(round(h[-1],4))
        P, T = P_to(W, rho, A_prop, V_to)
        RPM_to = RPM(T)
        V_tip = V_r(RPM_to,15.5,1)
        CT_to = CT(T,rho,A_prop,V_tip)
        FM_to = FM(CT_to, Cd0)
        eta_to = 0.866*0.98*FM_to*0.95  #eta_motor*eta_ESC*eta_prop
        E = np.append(E, E[-1]+P*dt/eta_to+E_sensors(dt))
        h = np.append(h,h[-1]+V_to*dt)
        t = np.append(t,t[-1]+dt)
        P_arr = np.append(P_arr,P)
   
    return E, E[-1], t, t[-1]-t_begin, P_arr

def E_climb(W, Cl_max, Cd_climb, S, phi, A_prop, t=np.array([0]), E=np.array([0]), P_arr=np.array([0]), h_cruise=500, h_trans=20):
    '''
    Calculates energy for climb
    '''
    t_begin = t[-1]
    h = np.array([h_trans])
    dt = 0.01

    while h[-1]<h_cruise:
        _,_,rho = isa(round(h[-1],4))
        P, V_ver, V_climb, T = P_climb(W, Cl_max, Cd_climb, S, phi, rho)
        V_a = v_axial_propeller(V_climb,T,rho,A_prop)
        eta_prop = n_prop(V_a,V_climb)
        eta_climb = 0.866*0.98*eta_prop*0.95  #eta_motor*eta_ESC*
        E = np.append(E, E[-1]+P*dt/eta_climb+E_sensors(dt))
        h = np.append(h,h[-1]+V_ver*dt)
        t = np.append(t,t[-1]+dt)
        P_arr = np.append(P_arr,P)
        
    return E, E[-1], t, t[-1]-t_begin, P_arr
        
def E_cruise(W, S, Cl_cruise, LD, phi, A_prop, t=np.array([0]), E=np.array([0]), P_arr=np.array([0]), h_cruise=500, h_trans=20, r=75000):
    '''
    Calculates energy for cruise
    '''
    t_begin = t[-1]
    phi = phi * np.pi/180
    s_climb = (h_cruise-h_trans)/np.sin(phi)
    s_glide = LD*(h_cruise-h_trans)
    s = r-s_climb-s_glide
    
    x = np.array([0])
    dt = 0.1
    
    while x[-1]<s:
        _,_,rho = isa(round(h_cruise,4))
        V_c = V_cruise(W, Cl_cruise, S, rho)
        Cd_c = Cl_cruise/LD
        q_c = q(rho, V_c) 
        P, T = P_cruise(q_c, S, Cd_c, V_c)
        V_a = v_axial_propeller(V_c,T,rho,A_prop)
        eta_prop = n_prop(V_a,V_c)
        eta_cruise = 0.866*0.98*eta_prop*0.95 #eta_motor*eta_ESC*eta_prop
        E = np.append(E, E[-1]+P*dt/eta_cruise+E_sensors(dt))
        x = np.append(x,x[-1]+V_c*dt)
        t = np.append(t,t[-1]+dt)
        P_arr = np.append(P_arr,P)
        
    return E, E[-1], t, t[-1]-t_begin, P_arr, s_glide

def E_gliding(s_glide, V_glide, t=np.array([0]), E=np.array([0]), P_arr=np.array([0])):
    '''
    Calculates energy for gliding
    '''
    t_begin = t[-1]
    x = np.array([0])
    dt = 0.1
    
    while x[-1]<s_glide:
        
        x = np.append(x,x[-1]+V_glide*dt)
        E = np.append(E, E[-1]+E_sensors(dt))
        t = np.append(t,t[-1]+dt)
        P_arr = np.append(P_arr,25.25)
        
    return E, E[-1], t, t[-1]-t_begin, P_arr

def E_landing(W, A_prop, Cd0, V_des, t=np.array([0]), E=np.array([0]), P_arr=np.array([0]), h_landing=20):
    '''
    Calculates energy for landing
    '''
    t_begin = t[-1]
    h = np.array([h_landing])
    dt = 0.01
    V_des = 4

    while h[-1]>0:
        _,_,rho = isa(round(h[-1],4))
        P, T = P_landing(W, A_prop,V_des)
        RPM_la = RPM(T)
        V_tip = V_r(RPM_la,15.5,1)
        CT_la = CT(T,rho,A_prop,V_tip)
        FM_la = FM(CT_la, Cd0)
        eta_landing = 0.866*0.98*FM_la*0.95  #eta_motor*eta_ESC*eta_prop
        E = np.append(E, E[-1]+P*dt/eta_landing+E_sensors(dt))
        h = np.append(h,h[-1]-V_des*dt)
        t = np.append(t,t[-1]+dt)
        P_arr = np.append(P_arr,P)
   
    return E, E[-1], t, t[-1]-t_begin, P_arr
        
    
def E_sensors(t, P_sensors=25.25):
    '''
    Calculates energy for sensors
    '''
    E_s = t*P_sensors
    
    return E_s


def E_trip(W, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, Cd0, h_cruise, V_des, t=np.array([0]), E=np.array([0]), P_arr=np.array([0]), lab=np.array([]), trip='go'):
    '''
    Calculates total energy for round trip
    '''
    
    E, E_T,t,t_t, P_arr = E_to(W, A_prop, Cd0, t, E, P_arr, V_to=6, h_trans=20)
    lab = np.append(lab, np.full((1, len(E)-len(lab)),'takeoff'+trip))
    E, E_Cl,t,t_Cl, P_arr = E_climb(W, Cl_max, Cd_climb, S, phi, A_prop, t, E, P_arr, h_cruise, h_trans=20)
    lab = np.append(lab, np.full((1, len(E)-len(lab)),'climb'+trip))
    E, E_Cr,t,t_Cr, P_arr, s_glide = E_cruise(W, S, Cl_cruise, LD, phi, A_prop, t, E, P_arr, h_cruise, h_trans=20, r=75000)
    lab = np.append(lab, np.full((1, len(E)-len(lab)),'cruise'+trip))
    E, E_g,t,t_g, P_arr = E_gliding(s_glide, 17, t, E, P_arr)
    lab = np.append(lab, np.full((1, len(E)-len(lab)),'glide'+trip))
    E, E_L,t,t_L, P_arr = E_landing(W, A_prop, Cd0, V_des, t, E, P_arr, h_landing=20)
    lab = np.append(lab, np.full((1, len(E)-len(lab)),'landing'+trip))
    
    return E, t, P_arr, lab


def plot_mission(m_bat_cell, e_d, EOL_corr, E_arr, t_arr, P_arr, lab):

    #SOC PLOT
    SoC_BOL_arr = (m_bat_cell*e_d-E_arr)/(m_bat_cell*e_d)*100
    SoC_EOL_arr = (m_bat_cell*e_d-E_arr/EOL_corr)/(m_bat_cell*e_d)*100
    
    SoC_fig, (SoC_BOL, SoC_EOL) = plt.subplots(1,2)
    
    SoC_BOL.plot(t_arr,SoC_BOL_arr)
    SoC_BOL.set_ylim(0,100)
    SoC_BOL.set_ylabel('SOC [%]')
    SoC_BOL.set_xlabel('Time [s]')
    SoC_EOL.plot(t_arr,SoC_EOL_arr)
    SoC_EOL.set_ylim(0,100)
    SoC_EOL.set_ylabel('SOC [%]')
    SoC_EOL.set_xlabel('Time [s]')
    
    #ENERGY PLOT
    E_fig, E_vs_t = plt.subplots()
    
    E_vs_t.plot(t_arr,E_arr)
    E_vs_t.set_ylabel('Energy drained from battery [J]')
    E_vs_t.set_xlabel('Time [s]')
    
#    E_vs_t.plot(t_arr[lab=='takeoffgo'], E_arr[lab=='takeoffgo'])
#    E_vs_t.plot(t_arr[lab=='takeoffback'], E_arr[lab=='takeoffback'])
#    E_vs_t.plot([],[],label="takeoff")
#    E_vs_t.plot(t_arr[lab=='climbgo'], E_arr[lab=='climbgo'])
#    E_vs_t.plot(t_arr[lab=='climbback'], E_arr[lab=='climbback'])
#    E_vs_t.plot([],[],label="climb")
#    E_vs_t.plot(t_arr[lab=='cruisego'], E_arr[lab=='cruisego'], label= "cruise")
#    E_vs_t.plot(t_arr[lab=='glidego'], E_arr[lab=='glidego'], label="glide")
#    E_vs_t.plot(t_arr[lab=='landinggo'], E_arr[lab=='landinggo'], label='landing')
#    plt.legend(loc="upper left", fontsize=22)


        
        
        