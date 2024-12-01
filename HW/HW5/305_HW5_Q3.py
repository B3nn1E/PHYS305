# -*- coding: utf-8 -*-
"""
Problem 3: Simulate the motion of 500 Brownian particles of radius 10 nm that 
interact by a Lennard-Jones potential.The particles should be in a 2D box of 
length 1 μm. Use a cutoff on the force so that the maximum value of the force
does not move the particles more than about a tenth of a particle diameter. 
Determine if the values you found in Problem 1 do indeed lead to aggregation, 
and whether values outside that range disperse. 
Provide the code you used and a set of values that lead to aggregation 
and one that doesn’t. Include a plot showing aggregation and one showing
dispersion.
"""

import numpy as np
import matplotlib.pyplot as plt

# Define the Lennard-Jones Force
def LJ_Force(r, sig, eps):
    # Force = dV / dr
    mask = r != 0 # make sure not to divided by 0
    force = np.zeros_like(r) # set force into 0 if r = 0
    force[mask] = 48 * eps * ((sig / r[mask])**13 - 0.5 * (sig / r[mask])**7)
    return force 

# Simulate Brownie dynamics 
def BrownianParticles(N, R, L, TotalTime, dt, Skip, eps, sig):
    kBT = 0.004
    eta = 0.001  
    zeta = 6 * np.pi * eta * R  
    c = zeta * np.sqrt(kBT)  
   
    Steps = round(TotalTime / dt / Skip)
   
    # Initialize a list to store positions
    x = np.zeros((2, N, Steps))
    x[:, :, 0] = np.random.rand(2, N) * L
   
    # Time stepping
    for i in range(1, Steps):
        for j in range(Skip):
            xi = c * np.random.randn(2, N)  # Stochastic force
            F_list = np.zeros((2, N))  # Initialize force list
           
            dx = x[0, :, i - 1].reshape(N, 1) - x[0, :, i - 1].reshape(1, N)
            dy = x[1, :, i - 1].reshape(N, 1) - x[1, :, i - 1].reshape(1, N)
            # Apply periodic boundary conditions
            dx -= L * np.round(dx / L)
            dy -= L * np.round(dy / L)
            r = np.sqrt(dx**2 + dy**2)
                      
            mask = (r < 3 * sig) & (r != 0) # avoid to interact with self so no 0 r
            F = LJ_Force(r, sig, eps) # calculate force
            F_list[0] = np.sum(mask * F * (dx / r), axis=1)
            F_list[1] = np.sum(mask * F * (dy / r), axis=1)
                            
            # Update positions
            x[:, :, i] += (dt / zeta) * (xi + F_list)
           
            # Apply periodic boundary conditions
            # Use '%' to get the remainder when divide by the length of box
            # Pretty much teleport the particle and throw it back in
            x[0, :, i] %= L 
            x[1, :, i] %= L

    return x

# Parameters (all in nm)
N = 500  
R = 10
L = 1000  
TotalTime = 1.0  
dt = 0.01  
Skip = 10  

# Aggregation vs. Dispersion
# Higher epsilon for agg so the particle can stick together
# Same particle so same sigma

eps_agg = 0.012 # 3kBT like my answer for 1
sig_agg = 20 # equal to the particles' diameters

eps_dis = 0.0001
sig_dis = 20

pos_agg = BrownianParticles(N, R, L, TotalTime, dt, Skip, eps_agg, sig_agg)
pos_dis = BrownianParticles(N, R, L, TotalTime, dt, Skip, eps_dis, sig_dis)

# Plot 
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Aggregation
for i in range(0, pos_agg.shape[2], 10):  # Plot every 10th frame
    axs[0].scatter(pos_agg[0, :, i], pos_agg[1, :, i], s=1)
axs[0].set_title('Aggregation of Brownian Particles (Lennard-Jones)')
axs[0].set_xlim(0, L)
axs[0].set_ylim(0, L)
axs[0].set_xlabel('x (nm)')
axs[0].set_ylabel('y (nm)')

# Dispersion
for i in range(0, pos_dis.shape[2], 10):
    axs[1].scatter(pos_dis[0, :, i], pos_dis[1, :, i], s=1)
axs[1].set_title('Dispersion of Brownian Particles (Lennard-Jones)')
axs[1].set_xlim(0, L)
axs[1].set_ylim(0, L)
axs[1].set_xlabel('x (nm)')
axs[1].set_ylabel('y (nm)')

plt.tight_layout()
plt.show()
