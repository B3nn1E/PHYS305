# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:53:41 2024

@author: kietb
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.animation import FFMpegWriter
from scipy.sparse import coo_array, eye_array
from scipy.sparse.linalg import spsolve

def SolveBZEq(L, N, TotalTime, dt, Skip, a, Du, Dv):
    # solves the BZ equation on a square, periodic grid of
    # Length L using N nodes. The equations are solved for a
    # total time, TotalTime, using a timestep dt. Data is plotted
    # and stored every skip time steps. The coefficient for the nonlinear
    # term in the Brusselator reaction is a. 
    # Dx and Dy are diffusion coefficients. 
    # Solve Belousov–Zhabotinsky reaction
    # Du/dt = D_u * Nu^2u + 1 - 4u + au^2v
    # Dv/dt = D_v * Nu^2v + 3u - au^2v
    
    # Create grid
    x, dx = np.linspace(0, L, N, retstep=True)
    y, dy = np.linspace(0, L, N, retstep=True)
    X, Y = np.meshgrid(x, y)
    
    # Number of time steps
    Steps = round(TotalTime / dt / Skip)
    
    # Total number of grid points
    Unknowns = N**2
    
    # Create identity matrix
    Iden = eye_array(Unknowns)
    
    # Create the Link matrix to number the nodes
    Link = np.array([i for i in range(Unknowns)])
    Link = np.reshape(Link, [N, N])
    
    # Define Laplacian matrix
    Lap = np.zeros((Unknowns, Unknowns))
    Lap[Link, Link] += -2 / dx**2 - 2 / dy**2
    Lap[Link, np.roll(Link, 1, axis=1)] += 1 / dx**2  # Left neighbor
    Lap[Link, np.roll(Link, -1, axis=1)] += 1 / dx**2  # Right neighbor
    Lap[Link, np.roll(Link, 1, axis=0)] += 1 / dy**2  # Bottom neighbor
    Lap[Link, np.roll(Link, -1, axis=0)] += 1 / dy**2  # Top neighbor
    
    Lap = coo_array(Lap)
    
    # Define left-hand side arrays
    Mu = Iden - (Du * dt) * Lap
    Mv = Iden - (Dv * dt) * Lap
    
    # Define initial conditions
    u = np.random.rand(N, N)
    v = np.random.rand(N, N)
    
    # Setup movie stuff
    metadata = dict(title='BZ reaction', artist='Matplotlib')
    writer = FFMpegWriter(fps=15, metadata=metadata)
    
    fig, ax = plt.subplots()
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    pcol = ax.pcolormesh(X, Y, u, cmap='jet')
    plt.colorbar(pcol, ax=ax, label="u concentration")
    
    # Begin Time Stepping
    with writer.saving(fig, 'BZ_reaction.mp4', 100):
        for n in range(Steps):
            for j in range(Skip):
                URHS = u + dt * (1 - 4 * u + a * u**2 * v)
                VRHS = v + dt * (3 * u - a * u**2 * v)
                
                u = spsolve(Mu, URHS.ravel()).reshape((N, N))
                v = spsolve(Mv, VRHS.ravel()).reshape((N, N))
            
            # Update plot and save frame
            pcol.set_array(u.ravel())
            ax.set_title(f"Step {n + 1}/{Steps}")
            writer.grab_frame()
            plt.pause(0.01)  # Allows the plot to refresh during execution

    return u, v


# Call the function
SolveBZEq(100, 100, 1, 0.01, 10, 4, 0.1, 0.1)