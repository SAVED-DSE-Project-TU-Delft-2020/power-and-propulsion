# -*- coding: utf-8 -*-
"""
Created on Thu May  7 12:19:59 2020

@author: Raven

This program is used to generate the loading diagram, which can be constructed from several cases.
run applyformat() to apply the formatting of the diagram
Altitude calculations based on ISA

Definition of Cd = k1 * CL^2 + k2 * CL + cd0

"""

from isacalc import isa
import matplotlib.pyplot as plt
import scipy as np




'''
Concept C: NACA35018
Cl/Cd @ Cldes = 50.0
Cmac = -0.0241
Alpha cruise = -0.024
Cl_max @ Re_transition = 0.6
Alpha stall @ Re_transition = 1.9
'''

#___CASE 1 CONSTANT ALITITUDE/SPEED CRUISE

def Case1(V, h, cd0, alpha, WS_range=np.arange(5,300,1)):
    '''
    calculates the wing loading for constant altitude and speed cruise
    designing for cruise speed
    '''
    T_W = []
    W_S = []
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    for WS in WS_range:
        TW = (q * cd0) / (alpha * WS)
        T_W.append(1/(V*TW))
#        T_W.append(TW)
        W_S.append(WS)
    
    plt.plot(W_S, T_W, label='Case 1', linewidth=4)


#___CASE 2 CONSTANT SPEED CLIMB
        
def Case2(cd0, k1, k2, V, h, alpha, beta, WS_range=np.arange(5,300,1)):
    '''
    calculates the wing loading for constant speed climb
    designing for climb 
    '''
    T_W = []
    W_S = []
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    for WS in WS_range:
        TW =  beta / alpha * (k1 * beta / q * WS + k2 + cd0 / (beta / q * WS) + 1)   
        T_W.append(1/(V*TW))
#        T_W.append(TW)
        W_S.append(WS)
        
    plt.plot(W_S, T_W, label='Case 2', linewidth=4)
    
    
#___CASE 3 CONSTANT ALTITUDE/SPEED TURN
    
def Case3(cd0, k1, k2, V, Rc, h, alpha, beta, WS_range=np.arange(5,300,1), g0=9.81):
    '''
    calculates the wing loading for constant altitude/ speed turn
    designing for turn radius
    '''
    T_W = []
    W_S = []
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    n   = np.sqrt(1+((V**2)/(g0*Rc))**2)
    
    for WS in WS_range:
        TW = beta / alpha * (k1 * n**2 * beta / q * WS + k2 * n + cd0 / (beta / q * WS))
        T_W.append(1/(V*TW))
#        T_W.append(TW)
        W_S.append(WS)

    plt.plot(W_S, T_W, label='Case 3', linewidth=4)

    
#___CASE 4 HORIZONTAL ACCELERATION

def Case4(cd0, k1, k2, dVdt, V, h, alpha, beta, WS_range=np.arange(5,300,1), g0=9.81):
    '''
    calculates the wing loading for horizontal acceleration
    designing for dVdt (acceleration)
    '''
    T_W = []
    W_S = []
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    for WS in WS_range:
        TW = beta / alpha * (k1 * beta / q * WS + k2 + cd0 / (beta / q * WS) + 1 / g0 * dVdt) 
        T_W.append(1/(V*TW))
#        T_W.append(TW)
        W_S.append(WS)

    plt.plot(W_S, T_W, label='Case 4', linewidth=4)

#___CASE 8 SERVICE CEILING
  
def Case8(cd0, k1, k2, V, dhdt, CL, h, alpha, beta, WS_range=np.arange(5,300,1)):
    '''
    calculates the wing loading for service ceiling
    designing for maximum altitude
    '''
    T_W = []
    W_S = []
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    for WS in WS_range:
        TW = beta / alpha * (k1 * beta / q * WS + k2 + cd0 / (beta / q * WS) + 1 / V * dhdt)
        T_W.append(1/(V*TW))
#        T_W.append(TW)
        W_S.append(WS)

    plt.plot(W_S, T_W, label='Case 8', linewidth=4)
    

#___CASE 0 STALL SPEED
    
def Case0(rho, Vs, CLmax):
    '''
    calculates the wing loading for stall speed
    designing for minimum airspeed
    '''    
    WS = 1/2 * rho * Vs**2 * CLmax
    plt.axvline(x=WS, label='Case 0', linewidth=4)
    print(WS)

    
#___CASE 00 ACCELERATION CLIMB
        
def Case00(cd0, k1, k2, dVdt, V, h, alpha, beta, WS_range=np.arange(5,300,1), g0=9.81):
    '''
    calculates the wing loading for accelerating climb
    designing for climb 
    '''
    T_W = []
    W_S = []
    
    _,_,rho = isa(h)
    q   = 0.5 * rho * V**2
    
    n = 1 + dVdt / g0
    
    for WS in WS_range:
        TW =  q / alpha * (cd0 / WS + k1 * (n * beta / q)**2 * WS + k2 * n * beta / q)  + beta / alpha * (1 + 1 / g0 * dVdt) 
        T_W.append(1/(V*TW))
#        T_W.append(TW)
        W_S.append(WS)
        
    plt.plot(W_S, T_W, label='Case 00', linewidth=4)    
    
def applyformat():
#    plt.xlim(0,260)
#    plt.ylim(0,1)
    plt.rcParams['font.size'] = 18.0
    plt.rc('axes', titlesize=28)
    plt.legend()
    plt.xlabel("W/S [$N/m^2$]")
    plt.ylabel("W/P [$N/W$]")
    plt.grid(color='grey', linestyle='-', linewidth=1)




'''

resulting WS:
WS = 122   

Values used for estimation:
    
cd0      = 0.0184
k1       = 0.0245
k2       = 0
    
Case1(28,500,cd0,0.9,WS_range)
Case2(cd0,k1,k1,6,500,0.9,1,WS_range)
Case3(cd0,k1,k2,20,9.81,7,500,0.9,1,WS_range)
Case4(cd0,k1,k2,5,9.81,10,0,0.9,1,WS_range)
Case8(cd0,k1,k2,28,5,0.6,2000,0.8,1,WS_range)
Case0(1.16727,14,1.06)
Case00(cd0,k1,k2,3,6,500,9.81,0.9,1,WS_range)
'''
    

'''
Concept A and D: NACA35120
Cl/Cd @ Cldes = 52.4
Cmac = 0.0113
Alpha cruise = 0.011
Cl_max @ Re_transition =1.8
Alpha stall @ Re_transition = 2.0

cd0 = 0.006
k1 = 0.06645
k2 = 0

Case1(28,500,0.006,0.9)
Case2(0.006,0.06645,0,6,500,0.9,1)
Case3(0.006,0.06645,0,20,11,500,0.9,1)
Case4(0.006,0.06645,0,5.4,28,500,0.9,1)
Case8(0.006,0.06645,0,28,5,0.011,2000,0.8,1)
Case0(1.16727,14,1.8)
Case00(0.006,0.06645,0,3,6,500,0.9,1)
'''
    



    
