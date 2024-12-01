# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 15:37:34 2024

@author: kietb
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
def Laplace1(H,W,Nx,Ny,VL,VR,VT,VB):
    
# solve Laplace’s equation on a rectangular domain
# of height H and width W. The number of gridpoints
# in x and y are Nx and Ny, respectively. The value of
# the electric potential on the boundaries is VL (left),
# VR (right), VT(top), and VB (bottom)


    x, dx = np.linspace(0, W, Nx, retstep = 'True')
    y, dy = np.linspace(0, H, Ny, retstep = 'True')
    
    X, Y = np.meshgrid(x, y)
    
    # total number of gridpoints where the solution will be computed
    Unknowns = Nx * Ny
    
    # create the Link matrix to number the nodes
    Link = np.array([i for i in range(Unknowns)])
    Link = np.reshape(Link, [Ny, Nx])
    
    Lap = np.zeros((Unknowns, Unknowns))
    
    Lap[Link, Link] += -2/dx ** 2 - 2/dy ** 2
    np.roll(Link, 1, axis = 1) # array to shift for link, how for to shift for 1, direction for axis
    Lap[Link, np.roll(Link, 1, axis = 1)] += 1/dx ** 2 # Left neighbor
    Lap[Link, np.roll(Link, - 1, axis = 1)] += 1/dx ** 2 # Right neighbor
    Lap[Link, np.roll(Link, 1, axis = 0)] += 1/dx ** 2 # Upper neighbor
    Lap[Link, np.roll(Link, - 1, axis = 0)] += 1/dx ** 2 # Lower neighbor
    
    # set a Dirichlet boundary condition on the top rows
    Lap[Link[0,:],:] = 0
    Lap[Link[0,:], Link[0,:]] = 1
    
    # set a Dirichlet boundary condition on the bottom rows
    Lap[Link[Ny - 1,:],:] = 0
    Lap[Link[Ny - 1,:], Link[Ny - 1,:]] = 1
    
    # set a Dirichlet boundary condition on the left rows
    Lap[Link[:,0],:] = 0
    Lap[Link[:,0], Link[:,0]] = 1
    
    # set a Dirichlet boundary condition on the right rows
    Lap[Link[:,Nx - 1],:] = 0
    Lap[Link[:, Nx - 1], Link[:,Nx - 1]] = 1

    # initialize source vector
    b = np.zeros(Unknowns)
    
    # set boundary conditions on source vector
    b[Link[0, :]] = VT
    b[Link[Ny - 1, :]] = VB
    b[Link[:, 0]] = VL
    b[Link[:, Nx - 1]] = VR
    
    # solve the equation
    Phi = np.linalg.solve(Lap, b)
    
    # reshape so that Phi has the same shape as the grid
    Phi = np.reshape(Phi, [Ny, Nx])
    
    # plot the solution
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    ax.plot_surface(X,Y,Phi,cmap=cm.jet)
    
    plt.show()
  
Laplace1(4, 2, 60, 40, 0, 5, 2, -1)

    
