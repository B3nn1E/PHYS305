# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:16:19 2024

@author: kietb
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

# Main function to perform reaction, diffusion, and advection
def ReactDiffuseAdvect(TotalTime, dt, Skip):
    
    D = 0.5  # Diffusion coefficient
    v = 0.25  # Advection velocity
    L = 2  # Length of the spatial domain
    
    # Calculate the number of steps for animation
    Steps = round(TotalTime / dt / Skip)
    
    N = 100  # Number of spatial points
    # Set up spatial grid and spacing (dx)
    x, dx = np.linspace(0, L, N, retstep=True, dtype='float')
    # Time array for plotting purposes
    Time = [Skip * dt * i for i in range(Steps)]
    
    # Initialize concentration array (2D) and second derivative array (1D)
    C = np.zeros((N, Steps))
    dCd2 = np.zeros(N)
    
    # Set up metadata and writer for saving animation
    metadata = dict(title='ReactDiffuseAdvect', artist='Matplotlib')
    writer = FFMpegWriter(fps=15, metadata=metadata)

    # Set up the figure and plot settings
    fig1, ax1 = plt.subplots()
    data, = ax1.plot([], [], 'b')
    ax1.set_xlim(x[0], x[-1])
    ax1.set_ylim(0, 1.5)
    ax1.set_title("Reaction-Diffusion-Advection System")
    ax1.set_xlabel("Position")
    ax1.set_ylabel("Concentration")

    # Start writing the animation
    with writer.saving(fig1, "ReactDiffuseAdvect.mp4", dpi=100):
        for i in range(1, Steps):
            # Update concentration array for the current time step
            C[:, i] = C[:, i - 1]
            
            # Perform reaction-diffusion-advection over `Skip` internal steps
            for j in range(Skip):
                # Boundary conditions at x=0
                dCd2[0] = D * 2 * (C[1, i - 1] - C[0, i - 1]) / (dx ** 2) + (2 - C[0, i - 1]) / dx
                # Interior points
                dCd2[1:N - 1] = D * (C[2:N, i - 1] - 2 * C[1:N - 1, i - 1] + C[0:N - 2, i - 1]) / (dx ** 2)
                # Boundary condition at x=L
                dCd2[N - 1] = D * 2 * (C[N - 2, i - 1] - C[N - 1, i - 1]) / (dx ** 2)
            
            # Calculate advection term
            dC = np.zeros(N)
            if v >= 0:
                # Upwind scheme for positive advection
                dC[0] = ((C[1, i - 1] - C[0, i - 1]) / dx) + (2 - C[0, i - 1]) / dx
                dC[1:N] = v * (C[1:N, i - 1] - C[:N - 1, i - 1]) / dx
            else:
                # Downwind scheme for negative advection
                dC[:N - 1] = v * (C[1:N, i - 1] - C[:N - 1, i - 1]) / dx
                dC[N - 1] = (1 + 0.5 * C[N - 1, i - 1])

            # Update concentration with reaction, diffusion, and advection
            C[:, i] = C[:, i - 1] + dt * (dC + dCd2 - C[:, i - 1])

            # Update plot for current time step in the animation
            data.set_data(x, C[:, i])
            writer.grab_frame()

    # Create a colormap plot for concentration over space and time
    fig2, ax2 = plt.subplots()
    cmap = plt.get_cmap('jet')
    im = ax2.pcolormesh(Time, x, C, shading='auto', cmap=cmap)
    fig2.colorbar(im, ax=ax2)
    ax2.set_title("Concentration over Time")
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Position")
    
    plt.show()
    return C, Time

# Example usage
TotalTime = 10
dt = 0.01
Skip = 10
C, Time = ReactDiffuseAdvect(TotalTime, dt, Skip)
