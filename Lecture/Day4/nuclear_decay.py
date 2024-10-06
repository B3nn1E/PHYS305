# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:13:12 2024

@author: kietb
"""
import numpy as np
import matplotlib.pyplot as plt

# def R dx/dt = x-x^2

def NuclearDecay(N0, k, T, Num): # add function as variable
    N = np.zeros(Num) # N = np.zeros(len(time))
    time, dt = np.linspace(0 , T, Num, 'retstep', 'true', 'float')
    # time = np.arrange (0, T, dt )
    N[0] = N0
    
    for i in range(1, Num): 
        N[i] = (1 - k*dt) * N[i - 1]
        # N[i] = N[i-1] + dt * R(i-1) (define a function R)
        
    return N, time

# define initial con
N0 = 10
k = 0.5
T = 5
Num = 100

N, time = NuclearDecay(N0, k, T, Num)

# define dt as a list of values  (np.concatenate, np.linspace)
# graph on multiple dt
plt.plot(time, N)


