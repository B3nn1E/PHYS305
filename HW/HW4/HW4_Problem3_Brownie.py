# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 19:00:27 2024

@author: kietb
"""

"""
Problem 3: 
    3. Simulate the 2D motion of a 2 µm diameter, neutrally-buoyant Brownian
    particle moving through water at room temperature using Langevin dynamics
    with a time step of 10-9 s and for a total time of 0.001 s. 
    Use the midpoint rule for handling the integration and compute the 
    randomforcing at the half time step. Show that the MSD as a function of 
    time is initially quadratic. At what time does it become linear? 
    Does this agree with what we discussed in class?
"""

import numpy as np
from numpy.random import randn
import matplotlib.pyplot as plt

def Brownie_Langevin (R, T, dt):
    # Constants
    m = 4 * np.pi * R**3/3 # mass
    kBT = 4e-14 # thermal NRG
    eta = 0.001 # viscosity
    zeta = 6 * np.pi * eta * R # drag coeff
    c = zeta * np.sqrt(kBT / m) # mag of stochasitc force
    
    steps = round(T / dt)
    
    # Initialize
    x = np.zeros((steps, 2))
    v = np.zeros((steps, 2))
    
    # Time stepping
    for i in range(1, steps):
        xi = c * np.random.randn(2) # random force
        
        # Midpoint velocity
        v[i, :] = ((1 - zeta * dt / (2 * m)) * v[i - 1, :] + dt * xi / m) / (1 + zeta * dt / (2 * m))
        v_mid = 0.5 * (v[i, :] + v[i - 1, :])
        
        # Update
        x[i, :] = x[i - 1, :] + dt * v_mid
    
    # Initialize
    Time = np.linspace(0, T, steps)
    MSD = np.zeros(steps)
    
    # Calculate MSD
    for i in range(steps):
        MSD[i] = np.mean(x[0, :i+1]**2 + x[1, :i+1]**2)
  
        
    # Plot
    plt.figure(figsize = (10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(x[0, :], x[1, :], label = 'Brownian 2D')
    plt.title('2D Motion of a Neutrally-Buoyant Brownian Particle')
    plt.grid(True)
    
    plt.subplot(1, 2, 2) 
    plt.plot(Time, MSD, label = 'Mean Squared Displacement')
    plt.xlabel = ('Time (s)')
    plt.ylabel('MSD (µm²)')
    plt.title('MSD vs Time')
    plt.grid(True)
    
    plt.legend()
    plt.show()

    return x, v, Time, MSD

# Call function
R = 1e-6
T = 0.001
dt = 1e-9

x, v, Time, MSD = Brownie_Langevin(R, T, dt)
    
    
        