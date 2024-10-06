# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 16:14:29 2024

@author: kietb
"""
import numpy as np
from numpy.random import randn
import matplotlib.pyplot as plt

def Langevin1(R, TotalTime, dt):
# simulates the dynamics of a neutrally-buoyant
# Brownian particle of raidus R moving thru water at room temp

    # parameters
    m = 4*np.pi*R**3/3 # mass of the particle
    kBT = 4*10**(-14) # thermal energy
    eta = 0.01 # viscosity of water
    zeta = 6*np.pi*eta*R # drag coeff
    c = zeta*(kBT/m)**(1/2) # magnitude of stochasitc force
    
    Steps = round(TotalTime/dt)
    
    #initialize lists for coordinates and velocities
    x = np.zeros((2, Steps))
    v = np.zeros((2, Steps))
    
    Time = [dt*i for i in range(Steps)]
    
    # begin time stepping
    for i in range(1, Steps):
        xi = c * np.random.rand(2)
        v[:,i] = ((1 - zeta*dt/2/m) * v[:,i-1] + dt*xi/m)/(1 + zeta*dt/2/m)
        vmid = 0.5 * (v[:,i] + v[:,i-1])
        x[:,i] = x[:,i-1] + dt*vmid
        
    plt.plot(x[0,:], x[1,:])
    
    return Time, x, v

Langevin1(10, 100, 1)
