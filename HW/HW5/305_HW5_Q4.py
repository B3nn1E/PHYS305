# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 23:25:08 2024

@author: kietb
"""

'''
Problem 4:
    Fisher’s equation: du/dt = D * d^2u/dx^2 + u * (1 - u)
describes a simple model for how a replicating species u spreads through 
space. This equation naturally leads to traveling wave solutions, 
where a sharp transition in u propagates at constant speed. 
Solve this equation on a domain of size 10, with boundary conditions 
u(x=0) = 1, and u(x=10) = 0 and an initial condition 
u(x,0) = 0.5(1 – tanh(x/0.1)). Determine how the propagation speed depends 
on D. You will need to determine how long to run the code in order to 
determine the propagation velocity. Provide a plot that shows the speed 
of the wave as a function of D.
'''

import numpy as np
import matplotlib.pyplot as plt

def fisher_double_deriv(dx, u):
    # Second derivative
    d2u_dx2 = np.zeros(len(u))
    for i in range(1, len(u) - 1):
        d2u_dx2[i] = (u[i + 1] - 2 * u[i] + u[i - 1]) / (dx * dx)
    return d2u_dx2


def fisher_main(dx, L, dt, TotalTime, D):
    N_steps = int(TotalTime / dt)
    N = int(L / dx) + 1
    x = np.linspace(0, L, N)
    u = np.zeros((N, N_steps))
    
    # Initial condition
    u[:, 0] = 0.5 * (1 - np.tanh((x - L / 2) / 0.1))
    
    # To store wave positions
    wave_positions = []

    for step in range(N_steps):
        d2u_dx2 = fisher_double_deriv(dx, u[:, step - 1])
        u[1:-1, step] = u[1:-1, step - 1] + dt * (D * d2u_dx2[1:-1] + u[1:-1, step - 1] * (1 - u[1:-1, step - 1]))
        
        # Boundary conditions
        u[0, step] = 1
        u[-1, step] = 0

        # Track wave front position
        wave_front = np.argmax(u[:, step] > 0.5) # return first position cross 0.5
        wave_positions.append(x[wave_front])

    speed = (wave_positions[-1] - wave_positions[0]) / TotalTime
    return speed, x, u

# Parameters
L = 10
dx = 0.1
TotalTime = 100
dt = 0.01

D_values = [0.1, 0.5, 1.0, 10.0, 20.0, 50.0, 100.0]
speeds = []

for D in D_values:
    speed, x, u = fisher_main(dx, L, dt, TotalTime, D)
    speeds.append(speed)

# Plot 
plt.figure()
plt.plot(D_values, speeds, marker='o')
plt.title("Wave Speed vs Diffusion Coefficient")
plt.xlabel("Diffusion Coefficient (D)")
plt.ylabel("Wave Speed")
plt.xscale('log')
plt.grid()
plt.show()
