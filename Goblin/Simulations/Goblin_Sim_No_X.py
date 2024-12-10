# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 23:38:02 2024

@author: kietb
"""
import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 4 * np.pi**2  # Gravitational constant (AU^3 / yr^2 / Msun)

# Masses (Msun)
M_sun = 1.0  # Sun
M_goblin = 1e-8  # Goblin

# Masses of outer planets (and Pluto)
masses = {
    "Uranus": 0.000044,  # Uranus
    "Neptune": 0.000052,  # Neptune
    "Pluto": 1.303e-4 * 3e-6,  # Pluto (~0.0022 Earth masses)
}

# Semi-major axes (AU)
a_goblin = 65  # Approximation for Goblin
ecc_goblin = 0.9  # Observed eccentricity for Goblin

# Semi major axes of outer planets (and Pluto)
semi_major_axes = {
    "Uranus": 19.2,
    "Neptune": 30.1,
    "Pluto": 39.5,
}

# Initial position (AU) at periapsis
r_goblin = np.array([a_goblin * (1 - ecc_goblin), 0, 0], dtype='float64')

# Initial velocity (AU/yr) using vis-viva equation
v_goblin = np.array([0, np.sqrt(G * M_sun * (2 / np.linalg.norm(r_goblin) - 1 / a_goblin)), 0], dtype='float64')

# Initialize planetary positions
planet_positions = {}
for planet, a in semi_major_axes.items():
    planet_positions[planet] = np.array([a, 0, 0], dtype='float64')

# Time parameters
T = 5200  # Time in years for the simulation
dt = 0.01  # Time step in years
steps = int(T / dt)

# Initialize arrays (position and velocity)
pos_goblin = np.zeros((steps, 3))
vel_goblin = np.zeros((steps, 3))
pos_goblin[0], vel_goblin[0] = r_goblin, v_goblin

# Gravitational force function
def grav_force(m1, r1, r2):
    r = np.linalg.norm(r1 - r2)
    if r == 0:
        return np.zeros(3)
    return G * m1 / r**3 * (r2 - r1)

# Simulation loop
for i in range(1, steps):
    # Force from the Sun
    F_goblin_sun = grav_force(M_sun, pos_goblin[i-1], np.zeros(3))

    # Add gravitational forces from outer planets
    F_goblin_outer = np.zeros(3)
    for planet, mass in masses.items():
        F_goblin_outer += grav_force(mass, pos_goblin[i-1], planet_positions[planet])

    # Net force
    F_goblin = F_goblin_sun + F_goblin_outer

    # Update velocity
    vel_goblin[i] = vel_goblin[i-1] + F_goblin * dt

    # Update position
    pos_goblin[i] = pos_goblin[i-1] + vel_goblin[i] * dt

# Plot Sun, the Goblin, and outer planets
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot([0], [0], [0], 'y*', label='Sun')
ax.plot(pos_goblin[:, 0], pos_goblin[:, 1], pos_goblin[:, 2], 'g-', label='The Goblin (2015 TG387)')

# Plot outer planets
planet_colors = {
    "Neptune": 'blue',
    "Uranus": 'teal',
    "Pluto": 'gray'
}
for planet, position in planet_positions.items():
    ax.scatter(position[0], position[1], position[2], color=planet_colors[planet], label=planet)

ax.set_title('Goblin Simulation with Outer Planets (No Planet X)')
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
ax.legend()
plt.show()