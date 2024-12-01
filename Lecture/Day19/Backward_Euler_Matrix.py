# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 16:22:11 2024

@author: kietb
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from Gauss_Jordan2 import GJ2

def DiffusionBE(D, L,N,TotalTime,dt,Skip):
    # solves the diffusion equation using Backward Euler
    
    # define parameters
    Nsteps = math.ceil(TotalTime/dt/Skip)
    
    # define grid
    x,dx = np.linspace(-L/2,L/2,N,retstep='True')
    
    Time = np.linspace(0,TotalTime,Nsteps)
    
    # initialize array for solution x
    C = np.zeros((N,Nsteps))
        
    # make BE matrix and identity matrix
    Lap = np.zeros((N,N))
    Iden = np.diag(np.ones(N),k=0)
    
    # use no flux boundary condition at x = 0 to set second derivative
    Lap[0,0] = -2
    Lap[0,1] = 2
    
    # define second derivative at interior points
    for i in range(1,N-1):
        Lap[i,i-1:i+2] = [1,-2,1]
        
        # use no flux boundary condition at x = L to set second derivative
        Lap[N-1,N-2] = 2
        Lap[N-1,N-1] = -2
        
        # define Lap to be the BE matrix
        Lap = Iden - (dt*D/dx**2)*Lap
        
    # setup movie plotting stuff
    fig1 = plt.figure()
    l, = plt.plot([],[],'b-')
    plt.xlim(0,L)
    plt.ylim(0,1)
    
    # define initial condition
    for j in range(N):
        C[j,0] = math.exp(-x[j]**2/0.1**2)
        
    # begin time stepping
    for i in range (1, Nsteps):
        C[:,i] = C[:, i - 1]
        for j in range(Skip):
            C[:, i] = GJ2(Lap, C[:, i])

        # plot it 
        l.set_data(x, C[:, i])
        plt.pause(0.02)
        
    cmap = plt.get_cmap('jet')
    
    fig2, ax2 = plt.subplots()
    im = ax2.pcolormesh(Time, x, C, shading = 'flat', cmap = cmap)
    
    plt.show()
    
    return Time, C
    
    