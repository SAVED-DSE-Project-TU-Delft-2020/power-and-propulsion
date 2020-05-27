# -*- coding: utf-8 -*-
"""
Created on Tue May 26 17:01:49 2020

@author: Raven
"""

from functions import E_trip, W_tot

A_prop = 0.4247
Cl_max = 1.539
Cd_climb = 0.088
S = 1.3
phi = 45
Cl_cruise = 0.2675095241
LD = 23.28
A = 6.923076923
e = 0.982
Cd0 = 0.0079
m_bat = 3.25
m_eng = 4.09
m_struc = 8.71
m_sensors = 0.363
eta = 0.73
DoD = 0.75
 
def main(m_bat, m_eng, m_struc, m_sensors, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, A, e, Cd0, eta, DoD):
    
    W_incl = W_tot(m_bat, m_eng, m_struc, m_sensors) 
    W_excl = W_tot(m_bat, m_eng, m_struc, m_sensors, PL=False)
    
    E_go, t_go = E_trip(W_incl, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, A, e, Cd0)
    E_back, t_back = E_trip(W_excl, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, A, e, Cd0)
    
    E_tot = E_go + E_back
    
    E_tot = E_tot/eta/DoD
    e = 180*3600
    
    mission_failpass = E_tot/e < m_bat
    
    return E_tot, mission_failpass