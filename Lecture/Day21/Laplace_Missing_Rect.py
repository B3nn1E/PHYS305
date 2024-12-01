# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 15:40:11 2024

@author: kietb
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
def LaplaceMissingRect(Hout,Wout,Hin,Win,Xin,Yin,Nx,Ny,V0):
# solve Laplace’s equation in a rectangular domain
# of height Hout and width Wout that has a rectangular
# conductor of size Hin by Win inside it with bottom left
# corner at (Xin,Yin). The number of gridpoints
# in x and y are Nx and Ny, respectively. The outer conductor
# is grounded and the inner one is at potential V0.

    # make grid
    x,dx = np.linspace(0,Wout,Nx,retstep='True')
    y,dy = np.linspace(0,Hout,Ny,retstep='True')
    X,Y = np.meshgrid(x,y)
    # mask out location of inner conductor
    InMask = (X >= Xin) & (X <= Xin + Win) & (Y >= Yin) & (Y <= Yin + Hin)
    
    # mask out where we do want to solve
    Mask = ~ InMask
    
    # find Edge of Inner Mask
    Edge = InMask & np.roll(Mask, 1, axis = 1)
    Edge = Edge | InMask & np.roll (Mask, -1, axis = 1)
    Edge = Edge | InMask & np.roll (Mask, 1, axis = 0)
    Edge = Edge | InMask & np.roll (Mask, -1, axis = 0)
    
    Mask = Mask | Edge
    Unknowns = np.count_nonzero(Mask)
    
    Link = -np.ones((Ny, Nx), dtype=int)
    Link[Mask] = np.array([i for i in range(Unknowns)], dtype = int)
    
    # define Laplacian matrix
    Lap = np.zeros((Unknowns, Unknowns))
    
    Del = [dy, dx] 
    
    for i in range(2):
        
        Here = Mask & np.roll(Mask, 1, axis = i)
        SLink = np.roll(Link, 1, axis = i)
        
        Lap[Link[Here], Link[Here]] += -1/Del[i]**2
        Lap[Link[Here], SLink[Here]] += 1/Del[i]**2
        Lap[SLink[Here], Link[Here]] += 1/Del[i]**2
        Lap[SLink[Here], SLink[Here]] += -1/Del[i]**2
        
    # define boundary conditions for top rows in Lap
    Lap[Link[0,:],:] = 0 
    Lap[Link[0,:], Link[0,:]] = 1
    
    # define boundary conditions for bottom rows in Lap
    Lap[Link[Ny - 1,:],:] = 0 
    Lap[Link[Ny - 1,:], Link[Ny-1,:]] = 1
    
    # define boundary conditions for left rows in Lap
    Lap[Link[:,0],:] = 0 
    Lap[Link[:,0], Link[:,0]] = 1
    
    # define boundary conditions for right rows in Lap
    Lap[Link[:,Nx - 1],:] = 0 
    Lap[Link[:,Nx - 1], Link[:,Nx - 1]] = 1
    
    # define boundary conditions for Edge rows in Lap
    Lap[Link[Edge],:] = 0 
    Lap[Link[Edge], Link[Edge]] = 1
    
    # initialize source vector
    b = np.zeros(Unknowns)
    
    # set boundary condition on the inner conductor
    b[Link[Edge]] = V0
    
    # initialize Phi matrix
    Phi = np.zeros((Ny, Nx))
    
    # solve the equation
    Phi[Mask] = np.linalg.solve(Lap, b)
    
    Phi[~Mask] = np.nan
    
    # Plot the solution
    #   light = Lightsource(90, 45)
    #   illuminated_surface = Light.shade(Phi, cmap = cm.jet, vmin = 0,, vmax = V0)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    ax.plot_surface(X, Y, Phi, antialiased = False, cmap= cm.jet, vmin = 0, vmax = V0)
    ax.set_xlim(-0.2, Wout + 0.2)
    ax.set_ylim(-0.2, Hout + 0.2)
    ax.set_zlim(-0.2, 1.2 * V0)
    
    plt.show()
    
        
LaplaceMissingRect(40, 10, 10, 4, 5, 2, 160, 60, 1)
    