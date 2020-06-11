# -*- coding: utf-8 -*-
"""
Created on Tue May 26 17:01:49 2020

@author: Raven
"""

from functions import E_trip, W_tot, plot_mission

#SAVED parameters

A_prop = 0.12173647185*4 #P&P propeller 15.5 inch 
Cl_max = 1.1371 #Aero
Cd_climb = 0.08113056706378069 #Aero
S = 1.318374324 #Planform
phi = 17 #ROC 5 m/s [deg]
Cl_cruise = 0.42059741185947325 #Aero
LD = 21.537381353841628 #Aero
Cd0 = 0.00976 #Aero
m_bat_cell = 18*0.175 #(P&P)
m_bat_casing = 0.45 #(P&P)
m_tot = 17.4754 #Interfacing
V_des = 4 #Turkish paper
DoD = 0.90 #Battery type
eta_bat = 0.99 #battery efficiency
EOL_corr = 0.75 #Kokam technical data sheet >3000 cycles
h_cruise = 500 #Mission profile


#Wingtra validation parameters
'''
A_prop = 0.07068583470577035*2 #Assumed based on figure 0.3 diameter
Cl_max = 1.1371 #Assumed equal
Cd_climb = 0.08113056706378069 #Assumed equal
S = 1.25*0.4 #Span * chord
phi = 17 #Assumed equal
Cl_cruise = 0.472887565 #Based on L=W 16 m/s cruise
LD = 21.537381353841628 #Assumed equal
Cd0 = 0.00976 #Assumed equal
V_des = 3 #Wingtra data
m_bat_cell = 0.604*0.875 #(P&P)
m_bat_casing = 0.604*0.125 #(P&P)
m_tot = 4.5 #Interfacing
DoD = 0.90 #Battery type
eta_bat = 0.99 #battery efficiency
EOL_corr = 0.75 #Kokam technical data sheet >3000 cycles
h_cruise = 500 #Mission profile
'''

def main(m_tot, m_bat_cell, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, Cd0, h_cruise, V_des, DoD, EOL_corr, eta_bat):
    '''
    Calculates the energy required based on aerodynamic parameters and the total mass.
    Calculates the actual battery mass required for the mission, and therefore returns a fail/pass for the total mission
    '''
    W_incl, W_excl = m_tot*9.81, (m_tot-3)*9.81
    
    E_arr, t_arr, P_arr, lab = E_trip(W_incl, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, Cd0, h_cruise, V_des)
    E_arr, t_arr, P_arr, lab = E_trip(W_excl, A_prop, Cl_max, Cd_climb, S, phi, Cl_cruise, LD, Cd0, h_cruise, V_des, t_arr, E_arr, P_arr, lab, trip='back')
    
    E_arr = E_arr/eta_bat
    
    E_tot_req = E_arr[-1]/DoD
    
    e_d = 246*3600 # for 3X6 KOKAM
#    e_d = 248*3600 # for 2X6 KOKAM
    
    plot_mission(m_bat_cell, e_d, EOL_corr, E_arr, t_arr, P_arr, lab)
    
    print(E_tot_req/e_d/EOL_corr +m_bat_casing , m_bat_cell+m_bat_casing)
    mission_failpass = E_tot_req/e_d/EOL_corr < m_bat_cell

    return E_tot_req, mission_failpass

