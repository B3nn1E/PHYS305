# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 15:35:25 2024

@author: kietb
"""

import numpy as np
import matplotlib.pyplot as plt 

def B(x,y):
    return np.array([0,0, np.sqrt(x**2 + y**2), 'float'])

def E(x,y):
    print(x, y)
    #return 0.01 * np.array([x, y, 0], 'float') / np.sqrt(x * x + y * y)
    return 0.01 * np.array([x, y, 0], 'float') / np.sqrt(x * x + y * y)

def phi(x,y):
    return 0.01 *1 / np.sqrt(x**2 + y**2)

r = np.array([0.9, 0, 0], 'float')
u = np.array([0.1, 0, 0], 'float')


dt = np.pi / 10
t_f = 500
steps = round(t_f / dt)

t = [i * dt for i in range(steps)]

m ,q, c = 1, 1, 1
traj = np.zeros((3, steps), 'float')

for i in range(steps):
    r_euler = r + dt * u
    E_mid = 0.5 * (E(r[0], r[1]) + E(r_euler[0], r_euler[1]))
    
    epsilon = (q / (2*m)) * E_mid
    
    u_min = u + epsilon * dt
    
    gamma = np.sqrt(1 + np.linalg.norm(u))**2 / (c**2)
    
    theta = (q * dt) / (m * gamma) * np.linalg.norm(B(r[0], r[1]))
    b_norm = B(r[0], r[1]) / np.linalg.norm(B(r[0], r[1]))
    
    u_min_par = np.dot(u_min, b_norm) * b_norm
    
    u_plu = u_min_par + (u_min - u_min_par) * np.cos(theta) + (np.cross(u_min, b_norm)) * np.sin(theta)
    
    u = u_plu + epsilon + dt
    
    r = dt * (r[i] + u / gamma)
    traj[:, i + 1] = r


plt.plot(traj[0, :], traj[1, :])