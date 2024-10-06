# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:42:47 2024

@author: kietb
"""

"""
Problem 1
The force from a real spring is nonlinear with an approximate force given by
    F = -kx + q/r^3 * r
where k and q are constants, r is the vector displacement from the origin 
with magnitude r. Use that k/m = 1, where m is the mass of the object 
connected to the spring. Solve for the two dimensional orbits in order to 
determine how the orbits depend on q/m. Use initial conditions
such that when q = 0, the orbit is circular.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

# Set up force 
def spring_force(r, q_m): # k/m = 1 => F/m = a = -r + q/m / r^3
    r_mag = np.linalg.norm(r) # r = sqrt(x^2 + y^2)
    return -r / r_mag + q_m * r / r_mag**3

# Create simulation
def spring(rs, vs, q_m, spring_force, dt, steps, writer, fig, ax):
    path = [] # empty list to store the plot direction later
    with writer.saving(fig, "2d_orbit.mp4", 100):
        for i in range(1, steps):
            F = spring_force(rs[i-1], q_m)
            a = F
            vs[i] = vs[i-1] + a * dt
            rs[i] = rs[i-1] + vs[i] * dt

            path.append(rs[i])
     # Plot       
    ax.plot(*zip(*path), '-', label=f'q/m = {q_m}')
        
# Initial conditions
if __name__ == '__main__':
    dt = 0.01
    T = 100
    steps = int(T / dt)
    
    rs = np.zeros((steps, 2))
    vs = np.zeros((steps, 2))
    
    r0 = np.array([1.0, 0.0])
    v0 = np.array([0.0, 1.0])
    rs[0] = r0
    vs[0] = v0
    
    # Test different values for q/m to see how orbit depends on q/m
    test_qms = [0.0, 0.5, 1.5, 2.5]

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    metadata = dict(title='2D Orbit', artist='Matplotlib')
    writer = FFMpegWriter(fps=15, metadata=metadata)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    # Call the function
    for q_m in test_qms:
        spring(rs, vs, q_m, spring_force, dt, steps, writer, fig, ax)

    ax.legend()
    fig.savefig('spring.png')
    plt.show()