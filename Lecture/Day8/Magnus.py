# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 15:32:54 2024

@author: kietb
"""

# Pick problem 1:
import numpy as np
import matplotlib.pyplot as plt

# Define the function
def magnus_force(m, v, g, A, rho, omega, s_m):
    return -g - (rho/m) * A * v**2 + s_m * omega * v

def RK4(m, v, g, A, rho, omega, s_m):
    k1 = dt * magnus_force(m, v, g, A, rho, omega, s_m)
    k2 = dt * magnus_force((v + k1 / 2),  g, A, rho, omega, s_m, m )
    k3 = dt * magnus_force((v + k2 / 2), g, A, rho, omega, s_m, m)
    k4 = dt * magnus_force((v + k3), g, A, rho, omega, s_m, m)

    k_final = v + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
    
    return k_final

def trajectory(m, v, g, A, rho, omega, s_m, t0, tf, dt):
    t = np.arrange(t0, tf, dt)
    v = np.zeros(len(t), 3)
    
    for i in range(1, len(t)):
        v[i] = RK4(m, v[i-1], g, A, rho, omega, s_m )
    return t, v 

if __name__ == '__main__':
    t0, tf, dt = [0, 10, 0.01]
    g = [0, 0, 9.8]
    m = 0.5
    A = 0.1
    s_m = 1
    rho = 1.293
    v0 = [10, 0, 0 ]
    omega = [0, 0, 10] 
    
    t, v = trajectory (m, v0, g, A, rho, omega, s_m, t0, tf,dt)
    fig, axs = plt.subplot(3,1)
    axs[0].plot(t, v[:,0])
    axs[1].plot(t, v[:,1])
    axs[2].plot(t, v[:,2])
    plt.show()
