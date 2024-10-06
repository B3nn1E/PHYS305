# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 16:34:01 2024

@author: kietb
"""

from numpy.random import rand
import matplotlib.pyplot as plt
import numpy as np

def RandomWalk1D(b , N):
# Computes a one-dimensional random walk of step-size b # and N steps

# Initialize position
    x = [0 for i in range(N)]
    for i in range(1, N):
        Flip = rand(1)
        
        if Flip < 0.5:
            x[i] = x [i -1] - b
        else:
            x[i] = x [i - 1] + b
            
# Plot result
    #plt.plot(x)

# Analyze the random walk (average displacement, MSD)

    AvgDisp = np.zeros(N)
    MSD = np.zeros(N)
    
    for i in range(N):
        for j in range(N-i):
            AvgDisp[i] += (x[i+j] - x[j]) / (N-i)
            MSD[i] += (x[i+j] - x[j])**2 / (N-i)
        
    
    plt.plot(AvgDisp, label = 'Average Displacement')
    plt.show()
    plt.plot(MSD, label = 'Mean Square Displacement')
    plt.show()
    plt.legend()
    return x, AvgDisp, MSD
    
RandomWalk1D(1, 1000)