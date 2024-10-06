# -*- coding: utf-8 -*-
"""
Problem 2: 
    Solve the 4-body problem for the motion of the Sun, Earth, Jupiter and 
Moon. Make a plot and a video showing the trajectories. Simulate long enough
that Jupiter does at least one full orbit around the Sun.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

# Define constants 
G = 4 * np.pi**2 # ( AU^3 / (yr^2 * Msun))

# Masses (Msun)
M_sun = 1.0 
M_earth = 3.003e-6 
M_jup = 9.545e-4  
M_moon = 3.694e-8 

# Initial positions relative to the sun (AU) and velocities (AU/yr)
# v = sqrt(GMsun/r)  = sqrt (4pi/r) = 2pi * sqrt (1/r)
r_sun = np.array([0, 0], dtype='float64')
v_sun = np.array([0, 0], dtype='float64')

r_earth = np.array([1, 0], dtype='float64')
v_earth = np.array([0, 2 * np.pi], dtype='float64') 


r_jup = np.array([5.2, 0], dtype='float64')
v_jup = np.array([0, 2 * np.pi * np.sqrt(1 / 5.2)], dtype='float64')

# For the moon, the positions and velocities will be relative to Earth
r_moon =  np.array([0.00257, 0], dtype='float64')
v_moon =  np.array([0, 2 * np.pi * np.sqrt(1 / 0.00257)], dtype='float64')


# Time parameters (in years)
T = 12  # Simulation time: 12 years (Jupiter's orbital period)
dt = 0.0001  # Time step (in years)
steps = int(T / dt)

# Arrays to store positions and velocities
pos_sun = np.zeros((steps, 2))
pos_earth = np.zeros((steps, 2))
pos_jup = np.zeros((steps, 2))
pos_moon = np.zeros((steps, 2))

pos_sun[0] = r_sun
pos_earth[0] = r_earth
pos_jup[0] = r_jup
pos_moon[0] = r_moon

speed_sun = np.zeros((steps, 2))
speed_earth = np.zeros((steps, 2))
speed_jup = np.zeros((steps, 2))
speed_moon = np.zeros((steps, 2))

speed_sun[0] = v_sun
speed_earth[0] = v_earth
speed_jup[0] = v_jup
speed_moon[0] = v_moon


# Calculate the gravitational force between 2 bodies
def grav_force(m1, m2, r1, r2):
    
    r = np.linalg.norm(r1 - r2) # distance between 2 bodies
    # r1 = [x1, y1] , r2 = [x2 , y2]
    # r = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    if r == 0: # prevent divided by 0 for the initial condition
        return np.zeros(2)
    
    F = G * m1 * m2 / r**3 * (r2 - r1)
    # Newton's Law of Gravitation: F = (G * m1 * m2) / r^2 * (r2 - r1)/ r (vector)
    return F

# Set up plot
fig, ax = plt.subplots(1, 1, figsize=(6, 6))
metadata = dict(title='4-Body Problem', artist='Matplotlib')
writer = FFMpegWriter(fps=15, metadata=metadata)

ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)

l1, = ax.plot([], [], 'y*', markersize=10)  # Sun
l2, = ax.plot([], [], 'bo', markersize=5)   # Earth
l3, = ax.plot([], [], 'o', color='orange', markersize=7)   # Jupiter
l4, = ax.plot([], [], 'ko', markersize=3)    # Moon

trace1, = ax.plot([], [], 'y-')
trace2, = ax.plot([], [], 'b-', alpha=0.5, linewidth=0.5)
trace3, = ax.plot([], [], '-', color='orange', alpha=0.5)
trace4, = ax.plot([], [], 'k-', alpha=0.5, linewidth=0.5)

# Set up's over, now simulation time
def Four_Body_Sim():
    with writer.saving(fig, "4body_simulation.mp4", 100):
        for i in range(1, steps):
            # Calculate gravitational forces on each body
            
            # Sun's forces
            F_sun_earth = grav_force(M_sun, M_earth, pos_sun[i-1], pos_earth[i-1])
            F_sun_jup = grav_force(M_sun, M_jup, pos_sun[i-1], pos_jup[i-1])
            F_sun_moon = grav_force(M_sun, M_moon, pos_sun[i-1], pos_earth[i-1]+pos_moon[i-1])

            # Earth's forces
            F_earth_sun = -F_sun_earth  # Newton's 3rd law
            F_earth_moon = grav_force(M_earth, M_moon, pos_earth[i-1], pos_moon[i-1])
            F_earth_jup = grav_force(M_earth, M_jup, pos_jup[i-1], pos_jup[i-1])

            # Jupiter's forces
            F_jup_sun = -F_sun_jup  
            F_jup_earth = -F_earth_jup
            F_jup_moon = grav_force(M_jup, M_moon, pos_jup[i-1], pos_moon[i-1])

            # Moon's forces
            F_moon_sun = -F_sun_moon
            F_moon_earth = -F_earth_moon
            F_moon_jup = -F_jup_moon

            # Net force on each body
            F_sun = F_sun_earth + F_sun_jup + F_sun_moon
            F_earth = F_earth_sun + F_earth_jup + F_earth_moon
            F_jup = F_jup_sun + F_jup_earth
            F_moon = F_moon_sun + F_moon_earth + F_moon_jup

            # Update velocity: v = adt = (F/m) * dt
            speed_sun[i] = speed_sun[i-1] + (F_sun / M_sun) * dt
            speed_earth[i] = speed_earth[i-1] + (F_earth / M_earth) * dt
            speed_jup[i] = speed_jup[i-1] + (F_jup / M_jup) * dt
            speed_moon[i] = speed_moon[i-1] + (F_moon / M_moon) * dt

            # Update position: r = vdt
            pos_sun[i] = pos_sun[i-1] + speed_sun[i] * dt
            pos_earth[i] = pos_earth[i-1] + speed_earth[i] * dt
            pos_jup[i] = pos_jup[i-1] + speed_jup[i] * dt
            pos_moon[i] = pos_moon[i-1] + speed_moon[i] * dt
            
            # Update the plot every steps taken
            if i % 100 == 0:
                l1.set_data(pos_sun[i, 0], pos_sun[i, 1])
                l2.set_data(pos_earth[i, 0], pos_earth[i, 1])
                l3.set_data(pos_jup[i, 0], pos_jup[i, 1])
                l4.set_data( pos_earth[i,0] + pos_moon[i, 0],  pos_earth[i,1] + pos_moon[i, 1])
        
                trace1.set_data(pos_sun[:i, 0], pos_sun[:i, 1])
                trace2.set_data(pos_earth[:i, 0], pos_earth[:i, 1])
                trace3.set_data(pos_jup[:i, 0], pos_jup[:i, 1])
                trace4.set_data(pos_earth[:i, 0]+pos_moon[:i, 1], pos_earth[:i,1] + pos_moon[:i, 1])
        
                writer.grab_frame()
                plt.pause(0.01) 

# Run the simulation
Four_Body_Sim()
plt.show()
