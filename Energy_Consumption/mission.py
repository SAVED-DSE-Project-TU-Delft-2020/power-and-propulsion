# -*- coding: utf-8 -*-
"""
Created on Tue May 26 17:01:49 2020

@author: Raven
"""

from functions import E_trip, W_tot, plot_mission

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
m_struc = 8.71 #Interfacing
m_sensors = 0.363 #Interfacing
eta = 0.63 #Battery type
DoD = 0.90 #Battery type
EOL_corr = 0.8 #Kokam technical data sheet >3000 cycles
h_cruise = 500 #Mission profile
 
def main(m_bat_cell, m_bat_casing, m_eng_prop, m_struc, m_sensors, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, A, e, Cd0, h_cruise, eta, DoD, EOL_corr):
    '''
    Calculates the energy required based on aerodynamic parameters and the total mass.
    Calculates the actual battery mass required for the mission, and therefore returns a fail/pass for the total mission
    '''
    W_incl, W_excl = W_tot(m_bat_cell+m_bat_casing, m_eng_prop, m_struc, m_sensors), W_tot(m_bat_cell+m_bat_casing, m_eng_prop, m_struc, m_sensors, PL=False) 
    
    E_go, t_go, P_go = E_trip(W_incl, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, A, e, Cd0, h_cruise)
    E_arr, t_arr, P_arr = E_trip(W_excl, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, A, e, Cd0, h_cruise, t_go, E_go, P_go)
    
    E_arr = E_arr/eta
    E_tot_req = E_arr[-1]/DoD
    
    e_d = 246*3600 # for 3X6 KOKAM
#    e_d = 248*3600 # for 2X6 KOKAM
    
    plot_mission(m_bat_cell, e_d, EOL_corr, E_arr, t_arr, P_arr)
    
    print(E_tot_req/e_d/EOL_corr +m_bat_casing , m_bat_cell+m_bat_casing)   
    mission_failpass = E_tot_req/e_d/EOL_corr < m_bat_cell

    return E_tot_req, mission_failpass



