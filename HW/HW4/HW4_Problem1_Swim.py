# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
"""
Problem 1: 
    A self-propelled object of mass m moves with velocity v = v0d, where v0
is the constant “swimming” speed of the object and d is a unit vector that 
points along the object’s direction of motion. When two of these objects 
are near each other they interact in such a way that the velocities and 
orientations try to entrain with each other. The equations of motion for 
these objects is governed by 4 equations. Write a function that will solve 
these equations given values for the constants.
"""

import numpy as np
import matplotlib.pyplot as plt

def swim_motion(m, zeta, alpha, beta, v0, dt, T):
    # Initial Condition
    v1 = np.array([v0, 0], dtype=float)
    v2 = np.array([0, v0], dtype=float)
    d1 = np.array([1, 0], dtype=float)   
    d2 = np.array([0, 1], dtype=float)   
    r1 = np.array([0, 0], dtype=float)
    r2 = np.array([1, 0], dtype=float)
    
    steps = int(T/dt)
    
    # Initialize lists
    v1_list = np.zeros((steps, 2))
    v2_list = np.zeros((steps, 2))
    d1_list = np.zeros((steps, 2))
    d2_list = np.zeros((steps, 2))
    
    for i in range (steps):
        r = np.linalg.norm(r2 - r1)
        
        # Using midpoint methods
        # Acceleration 
        dv1_dt = (zeta * (v0 * d1 - v1) + (alpha / r) * (v2 - v1))/m
        dv2_dt = (zeta * (v0 * d2 - v2) + (alpha / r) * (v1 - v2))/m
        
        # Midpoints
        v1_mid = v1 + 0.5 * dv1_dt * dt
        v2_mid = v2 + 0.5 * dv2_dt * dt

        dv1_dt_mid = (zeta * (v0 * d1 - v1_mid) + (alpha / r) * (v2_mid - v1_mid))/m
        dv2_dt_mid = (zeta * (v0 * d2 - v2_mid) + (alpha / r) * (v1_mid - v2_mid))/m
        
        # Update Velocity
        v1 += dv1_dt_mid * dt
        v2 += dv2_dt_mid * dt
        
        # Direction 
        # (a x b) x c = (a . c)b - (a . b)c
        # So dd1/dt = beta / r^2 ((d1 . d1)d2 - (d1 . d2)d1)
        # d1 . d1  = 1 so cancel that
        # Use np.dot, be smarter
        # Stupid whole, you should just np.cross twice
        dd1_dt = (beta / r**2) * (d2 - np.dot(d1, d2) * d1)
        dd2_dt = (beta / r**2) * (d1 - np.dot(d2, d1) * d2)
        
        # Midpoints
        d1_mid = d1 + 0.5 * dd1_dt * dt
        d2_mid = d2 + 0.5 * dd2_dt * dt
        
        dd1_dt_mid = (beta / r**2) * (d2 - np.dot(d1_mid, d2_mid) * d1)
        dd2_dt_mid = (beta / r**2) * (d1 - np.dot(d2_mid, d1_mid) * d2)
        
        # Update and normalized distance
        d1 += dd1_dt_mid * dt
        d2 += dd2_dt_mid * dt
        
        d1 /= np.linalg.norm(d1)
        d2 /= np.linalg.norm(d2)
        
        # Store
        v1_list[i] = v1
        v2_list[i] = v2
        d1_list[i] = d1
        d2_list[i] = d2
        
        # Update position, rinse and repeat
        r1 += v1 * dt
        r2 += v2 * dt
        
    return v1_list, v2_list, d1_list, d2_list

# Initial conditions
m = 1.0  
zeta = 0.1  
alpha = 0.5  
beta = 0.1
v0 = 1.0  
dt = 0.01
T = 10
        

v1s, v2s, d1s, d2s = swim_motion(m, zeta, alpha, beta, v0, dt, T)

# Plot
plt.plot (v1s[:, 0], label = 'v1')
plt.plot (v2s[:, 0], label = 'v2')
plt.xlabel('Time')
plt.ylabel('Velocity')
plt.legend()
plt.show()

