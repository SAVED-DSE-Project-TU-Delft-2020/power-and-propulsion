# -*- coding: utf-8 -*-
"""
Created on Fri May  8 11:13:48 2020

@author: Raven
"""
import math

def isa(h):
    #Constants
    
    T0=288.15
    R=287.
    g0=9.80665
    p0=101325.
    rho0=1.225
    
    if 0<=h<=11000:
        
        #Calculations
        
        T=T0-0.0065*h
        p=p0*(T/T0)**(-g0/(-0.0065*R))
        rho=p/(R*T)
        Tc=T-273.15
        Ppercent=p/p0*100.
        Rhopercent=rho/rho0*100.
      
        
    elif 11000<h<=20000:
    
        #Calculations
            
        T=T0-0.0065*11000
        p=22626*math.exp(-g0/(R*T)*(h-11000))
        rho=p/(R*T)
        Ppercent=p/p0*100
        Rhopercent=rho/rho0*100
        Tc=T-273.15

    elif 20000<h<=32000:
        
        #Calculations
                
        T=T0-0.0065*11000+0.001*(h-20000)
        p=5472*(T/216.65)**(-g0/(0.001*R))
        rho=p/(R*T)
        Ppercent=p/p0*100
        Rhopercent=rho/rho0*100
        Tc=T-273.15
            
    elif 32000<h<=47000:
        
        #Calculations
                
        T=T0-0.0065*11000+0.001*12000+0.0028*(h-32000)
        p=867*(T/228.65)**(-g0/(0.0028*R))
        rho=p/(R*T)
        Ppercent=p/p0*100
        Rhopercent=rho/rho0*100
        Tc=T-273.15
        
    elif 47000<h<=51000:
        
        #Calculations
            
        T=T0-0.0065*11000+0.001*12000+0.0028*15000
        p=111*math.exp(-g0/(R*T)*(h-47000))
        rho=p/(R*T)
        Ppercent=p/p0*100
        Rhopercent=rho/rho0*100
        Tc=T-273.15
    
    elif 51000<h<=71000:
        
        #Calculations                        
                            
        T=T0-0.0065*11000+0.001*12000+0.0028*15000-0.0028*(h-51000)
        p=67*(T/270.65)**(-g0/(-0.0028*R))
        rho=p/(R*T)
        Ppercent=p/p0*100
        Rhopercent=rho/rho0*100
        Tc=T-273.15
        
    elif 71000<h<=84852:
        
        #Calculations                        
                            
        T=T0-0.0065*11000+0.001*12000+0.0028*15000-0.0028*20000-0.002*(h-71000)
        p=67*(T/214.65)**(-g0/(-0.0028*R))
        rho=p/(R*T)
        Ppercent=p/p0*100
        Rhopercent=rho/rho0*100
        Tc=T-273.15


    return p, T, rho
