# -*- coding: utf-8 -*-
"""
Created on Thu May  7 15:14:58 2020

@author: Raven
"""

from cases import Case1, Case2, Case3, Case4, Case8
import scipy as np

cd0      = 0.0184
V        = 28
k1       = 0.0245
k2       = 0
g0       = 9.81
WS_range = np.arange(1,5,0.001)


def loading_diagram(cd0, k1, k2,q, V, dhdt, CL,g0,Rc, alpha, beta, WS_range):
    
    Case1(q, cd0, alpha, WS_range)
    Case2(cd0, k1, k2, V, dhdt, q, alpha, beta, WS_range)
    Case3(cd0, k1, k2, V, g0, Rc, q, alpha, beta, WS_range)
    Case4()
    Case8(cd0, k1, k2, V, dhdt, CL, q, alpha, beta, WS_range)
    
    
    
    