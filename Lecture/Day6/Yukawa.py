# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:29:40 2024

@author: kietb
"""
import matplotlib.pyplot as plt
import numpy as np

def yumo (m1, m2, T, dt, skip):
# compute the motions of 2 particles in Yukawa potential
    g = 1 # coupling constant
    lam = 0.5 # decay length 
    mu = m1*m2/(m1+m2) # reduced mass
    
# define step size 
    steps = round(T/dt/skip)
    
# initialize postions and velocity
    r1 = np.zeros((steps, 2), 'float')
    r2 = np.zeros((steps, 2), 'float')
    
    v1 = np.zeros((steps, 2), 'float')
    v2 = np.zeros((steps, 2), 'float')
    
    time = np.zeros((steps, 1), 'float')

# initial conditions
    r1[0,0] = -m2/(m1 + m2)
    r2[0,0] = m1/(m1 + m2)
    v1[0,1] = -np.sqrt(3 * mu* np.exp(-0.5)/2)/m1
    v2[0,1] = np.sqrt(3 * mu* np.exp(-0.5)/2)/m2
    
    
# initialize midpoint variables, don't need to be tracked in time
    r1mid = np.zeros((1, 2), 'float')
    r2mid = np.zeros((1, 2), 'float')
    v1mid = np.zeros((1, 2), 'float')
    v2mid = np.zeros((1, 2), 'float')
    
    from matplotlib.animation import FFMpegWriter
    metadata = dict(title = 'YukawaOrbit', artist = 'Matplotlib')
    writer = FFMpegWriter(fps = 15, metadata = metadata)
    
    fig1 = plt.figure()
    l1, = plt.plot([], [], 'bo')
    l2, = plt.plot([],[], 'ro')
    
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    
    with writer.saving(fig1, 'YukawaOrbit.mp4', 100):
          
    # start time stepping
        for i in range (1, steps):
            
            r1[i, :] = r1[i-1, :]
            r2[i, :] = r2[i-1, :]
            v1[i, :] = v1[i-1, :]
            v2[i, :] = v2[i-1, :]
            
            time[i] = time[i - 1]
            for j in range(skip):
                # compute r1 and r2 at 1/2 time step
                r1mid[0 , :] = r1[i, :] + dt* v1[i, :]/2
                r2mid[0 , :] = r2[i, :] + dt* v2[i, :]/2
                
                # define the force at the half time step
                r = np.sqrt((r1mid[0,0] - r2mid[0,0])**2 + (r1mid[0,1] - r2mid[0,1])**2)
                Fmag = g**2 * np.exp(-lam * r) * (lam * r +1)/ r**3
                
                # time step the velocities a full time step by midpoint
                v1mid[0, :] = 0.5 * v1[i, :]
                v2mid[0, :] = 0.5 * v2[i, :]
                
                v1[i, :] = v1[ i, :] + dt * Fmag * (r2[i, :] - r1[i, :])/m1
                v2[i, :] = v2[ i, :] + dt * Fmag * (r1[i, :] - r2[i, :])/m2
                
                v1mid[0, :] += 0.5 * v1[i, :]
                v2mid[0, :] += 0.5 * v2[i, :]
                
                # time step the positions a full time step by midpoint
                
                r1[i, :] = r1[i, :] + dt * v1mid[0, :]
                r2[i, :] = r2[i, :] + dt * v2mid[0, :]
            
            l1.set_data(r1[i, 0], r1[i, 1])
            l2.set_data(r2[i, 0], r2[i, 1])
            writer.grab_frame()
            plt.pause(0.02)
                
    #plt(r1[: , 0], r1[:, 1], 'b', r2[: , 0], r2[:, 1], 'r')
        
    return r1, v1, r2, v2, time


m1 = 1
m2 = 1
dt = 0.01 
T = 100
skip = 1
yumo (m1, m2, T, dt, skip)


"""
from matplotlib.animation import FFMpegWriter
metadata = dict(title = 'YukawaOrbit', artist = 'Matplotlib')
writer = FFMpegWriter(fps = 15, metadata = metadata)

fig1 = plt.figure()
l1, = plt.plot([], [], 'bo')
l2, = plt.plot([],[], 'ro')

plt.xlim(-1,1)
plt.ylim(-1,1)

with writer.saving(fig1, 'YukawaOrbit.mp4', 100):
    
    for i in range(1, steps):
        r1[i, :] = r1[i-1, :]
        r2[i, :] = r2[i-1, :]
        v1[i, :] = v1[i-1, :]
        v2[i, :] = v2[i-1, :]
         
        time[i] = time[i - 1]
        for j in range(skip):
        # compute r1 and r2 at 1/2 time step
            r1mid[0 , :] = r1[i, :] + dt* v1[i, :]/2
            r2mid[0 , :] = r2[i, :] + dt* v2[i, :]/2
             
        # define the force at the half time step
            r = np.sqrt((r1mid[0,0] - r2mid[0,0])**2 + (r1mid[0,1] - r2mid[0,1])**2)
            Fmag = g**2 * np.exp(-lam * r) * (lam * r +1)/ r**3
             
        # time step the velocities a full time step by midpoint
            v1mid[0, :] = 0.5 * v1[i, :]
            v2mid[0, :] = 0.5 * v2[i, :]
             
            v1[i, :] = v1[ i, :] + dt * Fmag * (r2[i, :] - r1[i, :])/m1
            v2[i, :] = v2[ i, :] + dt * Fmag * (r1[i, :] - r2[i, :])/m2
             
            v1mid[0, :] += 0.5 * v1[i, :]
            v2mid[0, :] += 0.5 * v2[i, :]
             
        # time step the positions a full time step by midpoint
             
            r1[i, :] = r1[i, :] + dt * v1mid[0, :]
            r2[i, :] = r2[i, :] + dt * v2mid[0, :]
        
        l1.set_data(r1[i, 0], r1[i, 1])
        l2.set_data(r2[i, 0], r2[i, 1])
        writer.grab_frame()
        plt.pause(0.02)
"""