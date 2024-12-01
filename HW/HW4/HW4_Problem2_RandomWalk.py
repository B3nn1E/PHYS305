# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 23:49:00 2024

@author: kietb
"""

"""
Problem 2:
    Simulate a 1000 step random walk in two dimensions using a fixed step 
    size. What is the average displacement over the entire walk? What is 
    the squared displacement as a function of the number of steps? 
    Does your results agree with what we predicted in class?
"""

import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt

def RandomWalk (b, N):
    # Initialize positions
    x = np.zeros(N)
    y = np.zeros(N)
    
    # Random walk
    for i in range (1, N):
        walk = rand(1)
        angle = rand(1) * 2 * np.pi # random direction can head to
        if walk < 0.5: # step forwards
            x[i] = x[i - 1] + b * np.cos(angle)
            y[i] = y[i - 1] + b * np.sin(angle)
        else: # step backwards
            x[i] = x[i - 1] - b * np.cos(angle)
            y[i] = y[i - 1] - b * np.sin(angle)

    # MSD and AvgDisp
    MSD = np.zeros(N)
    AvgDisp = np.zeros(N)
    
    for i in range (N):
        for j in range (N-i):
            displacement = np.sqrt((x[i+j] - x[j])**2 + (y[i+j] - y[j])**2)
            MSD[i] += displacement**2 / (N-i)
            AvgDisp[i] += displacement / (N-i)
            
    # Plot
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(x, y, label = 'Random Walk Path')
    plt.xlabel('x position')
    plt.ylabel('y position')
    plt.title('Random Walk 2D')
    
    plt.subplot(1, 2, 2)
    plt.plot(MSD, label='Mean Squared Displacement')
    plt.plot(AvgDisp, label='Average Displacement')
    plt.xlabel('Number of Steps')
    plt.ylabel('Displacement')
    plt.title('Mean Squared Displacement vs Average Displacement')
    
    plt.legend()
    plt.show()
    
    return x, y, MSD, AvgDisp

# Initial Conditions
b = 1 # Step size
N = 1000 # Step #
x , y, MSD, AvgDisp = RandomWalk(b, N)