# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 18:47:51 2024

@author: kietb
"""

"""
Problem 1:
    
   Use the forward Euler method to solve: 
   dx/dt = x - x^2 with initial condition, x(0) = 0.1. Integrate 
from t = 0 to t = 30.

   Investigate how the solution depends on time step, using a range of values
between ∆t = 0.01 and 3. Pay special attention to the range of ∆t 
between 2 and 3. Provide graphs to show all the behaviors you find. 
Try to come up with an explanation for what happens between 
∆t = 2 and 3. 
"""
import numpy as np
import matplotlib.pyplot as plt

# define the function x - x^2
def funcR (x):
    return x - x**2

# apply Forward Euler method
def ForwardEuler (N0, T, dt):
    
    time = np.arange(0, T, dt)
    N = np.zeros(len(time))
    N[0] = N0
    
    for i in range(1, len(time)): 
        N[i] = N[i-1] + dt * funcR(N[i-1])
        
    return N, time

# define initial con and plot
N0 = 0.1
T = 30

dt_list = [0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

for dt in dt_list: 
    N, time = ForwardEuler (N0, T, dt)
    plt.plot(time, N, label = f'dt = {dt}')

plt.xlabel('Time (t)')
plt.ylabel('x')
plt.title('Solve the equation dx/dt = x - x^2 forward Euler method')
plt.legend()
plt.grid()
plt.show()


