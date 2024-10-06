# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 15:55:18 2024

@author: kietb
"""

import numpy.random as npr
import numpy as np
import matplotlib.pyplot as plt

def BrownianParticles(N,R,L,TotalTime,dt,Skip):
# simulates the Brownain dynamics of N neutrally-buoyant
# particles of radius R moving through water at room
# temperature and interacting via Lennard-Jones potentials.
# Periodic boundaries are used with a box size of L.

    # parameters
    kBT = 0.004 # thermal energy in pN*um
    eta = 0.001 # viscosity of water
    zeta = 6*np.pi*eta*R # drag coefficient
    c = zeta*(kBT)**(1/2) # magnitude of stochastic force
    
    Steps = round(TotalTime/dt/Skip)
    
    # initialize lists to store positions and time
    x = np.zeros((2, N, Steps))
    Time = [Skip*dt*i for i in range(Steps)]
    
    # initialize matrix to store the particle-particle interactions
    Fi = np.zeroes((2, N))
    
    for i in range (1, Steps):
        x[:,:,i] = x[:,:,i-1]
        
        for j in range (Skip):
            xi = c * np.random.randn(2, N)
            
            """
            
            for m in range(N):
                for n in range(N):
                    
                    Dx[m,n] = x[0,m,i] - x[0,n,i]
                    Dy[m,n] = x[1,m,i] – x[1,n,i]
            """ # shit way
            
            # Create matrix
            Dx = np.meshgrid(x[0,:,i-1])
            Dx = Dx - np.transpose(Dx)
            
            Dy = np.meshgrid(x[1,:,i])
            Dy = Dy - np.transpose(Dy)
            
            Dx[Dx<-L/2] = Dx[Dx<-L/2] + L
            Dx[Dx>L/2] = Dx[Dx>L/2] - L
            
            Dy[Dy<-L/2] = Dy[Dy<-L/2] + L
            Dy[Dy>L/2] = Dy[Dy>L/2] - L

            r = np.sqrt(Dx**2 + Dy**2)
            
            for j in range(N):
                r[j,j] = 10
            
            FiMag = np.e*( (sigma/r)**(14) - (sigma/r)**8 )
            FiMag*Dx (x-component)
            FiMag*Dy (y-component)

                    
            