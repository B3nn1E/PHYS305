# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 16:45:36 2024

@author: kietb
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def Laplace1(H, W, Nx, Ny, VL, VR, VT, VB):
    """
    Solve Laplace's equation in a rectangular domain.

    Parameters:
    H (float): Height of the domain
    W (float): Width of the domain
    Nx (int): Number of grid points in the x-direction
    Ny (int): Number of grid points in the y-direction
    VL, VR, VT, VB (float): Boundary potentials (left, right, top, bottom)

    Returns:
    phi (2D array): Computed electric potential on the grid
    """
    # Create grid
    x, dx = np.linspace(0, W, Nx, retstep=True)
    y, dy = np.linspace(0, H, Ny, retstep=True)
    X, Y = np.meshgrid(x, y)

    unknowns = Nx * Ny
    Link = np.arange(unknowns).reshape(Ny, Nx)

    # Define Laplacian matrix
    Lap = np.zeros((unknowns, unknowns))
    Lap[Link, Link] += -2 / dx**2 - 2 / dy**2
    Lap[Link, np.roll(Link, 1, axis=1)] += 1 / dx**2
    Lap[Link, np.roll(Link, -1, axis=1)] += 1 / dx**2
    Lap[Link, np.roll(Link, 1, axis=0)] += 1 / dy**2
    Lap[Link, np.roll(Link, -1, axis=0)] += 1 / dy**2

    # Apply boundary conditions
    Lap[Link[0, :], :] = 0
    Lap[Link[0, :], Link[0, :]] = 1
    Lap[Link[Ny-1, :], :] = 0
    Lap[Link[Ny-1, :], Link[Ny-1, :]] = 1
    Lap[Link[:, 0], :] = 0
    Lap[Link[:, 0], Link[:, 0]] = 1
    Lap[Link[:, Nx-1], :] = 0
    Lap[Link[:, Nx-1], Link[:, Nx-1]] = 1

    b = np.zeros(unknowns)
    b[Link[0, :]] = VT
    b[Link[Ny-1, :]] = VB
    b[Link[:, 0]] = VL
    b[Link[:, Nx-1]] = VR

    phi = np.linalg.solve(Lap, b).reshape(Ny, Nx)

    # Plot solution
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(X, Y, phi, cmap=cm.jet)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Electric Potential')
    plt.show()

    return phi


