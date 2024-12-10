# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 15:12:08 2024

@author: kietb
"""
'''
Problem 1: 
    A rectangular space defined by 0 < x < 4 and 0 < y < 2 has a 
    spatially non-uniform charge density given by 
    ρ(x,y) = sin(2πx)sin2(πy). 
    Determine the electric potential inside the domain assuming 
    that the surface of the rectangular region is grounded.
'''
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

def SolveLaplace(Lx, Ly, Nx, Ny, VL, VR, VT, VB):
    """
    Solve Laplace's equation for a rectangular domain.

    Parameters:
    Lx, Ly (float): Width and height of the domain
    Nx, Ny (int): Number of grid points in the x- and y-directions
    VL, VR, VT, VB (float): Boundary potentials (left, right, top, bottom)

    Returns:
    phi (2D array): Computed electric potential on the grid
    """
    # Create grid
    x, dx = np.linspace(0, Lx, Nx, retstep=True)
    y, dy = np.linspace(0, Ly, Ny, retstep=True)
    X, Y = np.meshgrid(x, y)

    unknowns = Nx * Ny
    Link = np.arange(unknowns).reshape(Ny, Nx)

    # Define Laplacian matrix
    Lap = np.zeros((unknowns, unknowns))
    Lap[Link, Link] += -4 / dx**2 - 4 / dy**2
    Lap[Link, np.roll(Link, 1, axis=1)] += 1 / dx**2  # Left neighbor
    Lap[Link, np.roll(Link, -1, axis=1)] += 1 / dx**2  # Right neighbor
    Lap[Link, np.roll(Link, 1, axis=0)] += 1 / dy**2  # Bottom neighbor
    Lap[Link, np.roll(Link, -1, axis=0)] += 1 / dy**2  # Top neighbor

    # Apply boundary conditions
    Lap[Link[0, :], :] = 0
    Lap[Link[0, :], Link[0, :]] = 1
    Lap[Link[Ny-1, :], :] = 0
    Lap[Link[Ny-1, :], Link[Ny-1, :]] = 1
    Lap[Link[:, 0], :] = 0
    Lap[Link[:, 0], Link[:, 0]] = 1
    Lap[Link[:, Nx-1], :] = 0
    Lap[Link[:, Nx-1], Link[:, Nx-1]] = 1

    # Define source term
    rho = np.sin(2 * np.pi * X) * (np.sin(np.pi * Y)**2)  # Charge density
    b = -rho.flatten() * dx * dy  # Source vector
    b[Link[0, :]] = VT  # Top boundary
    b[Link[Ny-1, :]] = VB  # Bottom boundary
    b[Link[:, 0]] = VL  # Left boundary
    b[Link[:, Nx-1]] = VR  # Right boundary

    # Solve for phi
    phi = np.linalg.solve(Lap, b).reshape(Ny, Nx)

    # Plot solution
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(X, Y, phi, cmap=cm.jet)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Electric Potential in Rectangular Space')
    plt.show()

    return phi

# Example usage
Lx, Ly = 4, 2  # Domain dimensions
Nx, Ny = 50, 25  # Grid points
VL, VR, VT, VB = 0, 0, 0, 0  # Boundary potentials

phi = SolveLaplace(Lx, Ly, Nx, Ny, VL, VR, VT, VB)
