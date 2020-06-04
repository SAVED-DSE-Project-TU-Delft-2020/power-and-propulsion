# -*- coding: utf-8 -*-
"""
Created on Tue May 26 17:01:49 2020

@author: Raven
"""

from functions import E_trip, W_tot

A_prop = 0.456 #P&P propeller
Cl_max = 1.539 #Aero
Cd_climb = 0.088 #Aero
S = 1.3 #Aero
phi = 17 #ROC 5 m/s
Cl_cruise = 0.2675095241 #Aero
LD = 23.28 #Aero
A = 6.923076923 #Aero
e = 0.982 #Aero
Cd0 = 0.0079 #Aero
m_bat_cell = 18*0.175 #Interfacing (P&P)
m_bat_casing = 0.45 #Interfacing (P&P)
m_eng_prop = 4*0.4+0.5 #Interfacing
m_struc = 8.71+0.5 #Interfacing
m_sensors = 0.363 #Interfacing
eta = 0.6 #Battery type
DoD = 0.90 #Battery type
EOLcorr = 0.75 #Kokam technical data sheet >3000 cycles
 
def main(m_bat_cell, m_bat_casing, m_eng_prop, m_struc, m_sensors, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, A, e, Cd0, eta, DoD):
    '''
    Calculates the energy required based on aerodynamic parameters and the total mass.
    Calculates the actual battery mass required for the mission, and therefore returns a fail/pass for the total mission
    '''
    W_incl, W_excl = W_tot(m_bat_cell+m_bat_casing, m_eng_prop, m_struc, m_sensors), W_tot(m_bat_cell+m_bat_casing, m_eng_prop, m_struc, m_sensors, PL=False) 
    
    E_go, t_go = E_trip(W_incl, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, A, e, Cd0)
    E_back, t_back = E_trip(W_excl, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, A, e, Cd0)
    
    E_tot = E_go + E_back
    
    E_tot = E_tot/eta/DoD
    e_d = 246*3600 # for 3X6 KOKAM
#    e_d = 248*3600 # for 2X6 KOKAM
    print(E_tot/e_d/EOLcorr +m_bat_casing , m_bat_cell+m_bat_casing)
    mission_failpass = E_tot/e_d/EOLcorr < m_bat_cell

    return E_tot, mission_failpass

