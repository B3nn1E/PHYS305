import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 4 * np.pi**2  # Gravitational constant (AU^3 / yr^2 / Msun)

# Masses (Msun)
M_sun = 1.0  # Sun
M_goblin = 1e-8  # Goblin
M_Earth = 3e-6 # Earth 
M_planetX = 10 * M_Earth # Planet X (~10 Earth masses)

# Masses of outer planets (and Pluto)
masses = {
    "Uranus": 0.000044, # Uranus
    "Neptune": 0.000052, # Neptune
    "Pluto": 1.303e-4 * M_Earth,  # Pluto (~0.0022 Earth masses)
}


# Semi-major axes (AU)

# From (Sheppard et al., 2018)
a_goblin = 65 # Approximation for Goblin
ecc_goblin = 0.9  # Observation eccentricity for Goblin

# From (Batygin and Brown, 2016)
a_planetX = 300  # Approximation for Planet X 
ecc_planetX = 0.6  # Speculated eccentricity for Planet X

# Semi major axes of outer planets (and Pluto)
semi_major_axes = {
    "Uranus": 19.2,
    "Neptune": 30.1,
    "Pluto": 39.5,
}

# Initial positions (AU) at periapsis
r_goblin = np.array([a_goblin * (1 - ecc_goblin), 0, 0], dtype='float64')
r_planetX = np.array([a_planetX * (1 - ecc_planetX), 0, 0], dtype='float64')

# Initial velocities (AU/yr) using vis-viva equation
v_goblin = np.array([0, np.sqrt(G * M_sun * (2 / np.linalg.norm(r_goblin) - 1 / a_goblin)), 0], dtype='float64')
v_planetX = np.array([0, np.sqrt(G * M_sun * (2 / np.linalg.norm(r_planetX) - 1 / a_planetX)), 0], dtype='float64')

# Initialize planetary positions
planet_positions = {}
for planet, a in semi_major_axes.items():
    planet_positions[planet] = np.array([a, 0, 0], dtype='float64')

# Time parameters
T = 5200  # Time in years for Planet X to complete 1 orbit
dt = 0.01  # Time step in years
steps = int(T / dt)

# Initialize arrays (position and velocity)
pos_goblin = np.zeros((steps, 3))
pos_planetX = np.zeros((steps, 3))
vel_goblin = np.zeros((steps, 3))
vel_planetX = np.zeros((steps, 3))
pos_goblin[0], vel_goblin[0] = r_goblin, v_goblin
pos_planetX[0], vel_planetX[0] = r_planetX, v_planetX

# Gravitational force function
def grav_force(m1, r1, r2):
    r = np.linalg.norm(r1 - r2)
    if r == 0:
        return np.zeros(3)
    return G * m1 / r**3 * (r2 - r1)

# Simulation loop
for i in range(1, steps):
    # Forces from Sun and Planet X
    F_goblin_sun = grav_force(M_sun, pos_goblin[i-1], np.zeros(3))
    F_goblin_planetX = grav_force(M_planetX, pos_goblin[i-1], pos_planetX[i-1])
    F_planetX_sun = grav_force(M_sun, pos_planetX[i-1], np.zeros(3))
    F_planetX_goblin = -F_goblin_planetX
    '''
    # Add gravitational forces from outer planets
    F_goblin_outer = np.zeros(3)
    for planet, mass in masses.items():
        F_goblin_outer += grav_force(mass, pos_goblin[i-1], planet_positions[planet])
    '''
    # Net forces
    F_goblin = F_goblin_sun + F_goblin_planetX 
    F_planetX = F_planetX_sun + F_planetX_goblin

    # Update velocities
    vel_goblin[i] = vel_goblin[i-1] + F_goblin * dt
    vel_planetX[i] = vel_planetX[i-1] + F_planetX * dt

    # Update positions
    pos_goblin[i] = pos_goblin[i-1] + vel_goblin[i] * dt
    pos_planetX[i] = pos_planetX[i-1] + vel_planetX[i] * dt

# Simulation time!

# Plot Sun, the Goblin and Planet X's orbits
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot([0], [0], [0], 'y*', label='Sun')
ax.plot(pos_goblin[:, 0], pos_goblin[:, 1], pos_goblin[:, 2], 'g-', label='The Goblin (2015 TG387)')
ax.plot(pos_planetX[:, 0], pos_planetX[:, 1], pos_planetX[:, 2], 'r-', label='Planet X')
'''
# Plot outer planets
planet_colors = {
    "Neptune": 'blue',
    "Uranus": 'teal',
    "Pluto": 'gray'
}
for planet, position in planet_positions.items():
    ax.scatter(position[0], position[1], position[2], color = planet_colors[planet], label=planet)
'''
ax.set_title('Goblin and Planet X Simulation (No Outer Planets)')
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
ax.legend()
plt.show()