# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 15:36:34 2024

@author: kietb
"""

# Advection Equation dC/dt = -v * dC/dx = - d(vC)/dx
# Continuity Equation dC/dt = -dJ/dx
# Periodic boundary condition
# Stablility condition d(x^2)> 2* D * dt

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

def AdvectionCD(v, L, N, TotalTime, dt):
    # simulates the motion of particles that are driven by a
    # velocity v and are contained in a region of size L
    # Periodic boundary condition are used 
    
    # define parameters
    Nsteps = round(TotalTime/dt) # Numbers of time steps
    
    # define grids
    x, dx = np.linspace(-L/2, L/2, retstep = 'True')   
    Time = np.linspace(0,TotalTime,Nsteps)
    
    # initialize concentration matrix
    C = np.zeros((N, Nsteps))
    
    # define initial condition
    for j in range (N):
        C[j, 0] = 0.5 * (np.tanh(x[j]/(0.1 * L)) + 1)
        
    # initialize first derivative array
    dC = np.zeros(N)
    
    # setup movies stuff
    metadata = dict(title = 'Avect', artist = 'Matplotlib')
    writer = FFMpegWriter(fps = 15, metadata = metadata)
    
    fig1 = plt.figure()
    l, = plt.plot([], [], 'b-')
    plt.xlim(-L/2, L/2)
    plt.ylim(0, 1)
    
    # begin time stepping
    with writer.saving(fig1, "AdvectionCD.mp4", 100):
        for i in range(1, Nsteps):
            if v >= 0:
            # compute first derivative
                dC[0] = (C[0, i - 1] - C[N - 1, i - 1]) /  dx
                dC[1:N] = (C[1:N, i - 1] - C[:N - 1, i - 1]) / dx
            else:
                dC[:N-1] = (C[1:N,i-1] - C[:N-1,i-1])/dx
                dC[N-1] = (C[0,i-1] - C[N-1,i-1])/dx
                
            # time step concentration
            C[:, i] = C[:, i - 1]  - dt * v * dC[:]
            
            l.set_data(x, C[:, i])
            writer.grab_frame()
            plt.pause(0.01)
            
        cmap = plt.get_cmap('inferno')
        
        fig2, ax2 = plt.subplots()
        im = ax2.pcolormesh(Time, x,C, cmap = cmap)
        
        plt.show()
        
        return Time, C
    
T ,C = AdvectionCD(1,2,50,1,0.001)

