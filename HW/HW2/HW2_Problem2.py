# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 16:34:24 2024

@author: kietb
"""

"""
Problem 2:
    ...If the concentrations of A and B are held fixed by external processes, the Law of Mass Action
leads to the following form for the dynamics of x and y:
    dx/dt = a-x + x^2*y
    dy/dt = b - x^2*y

where a and b are constants. Write a function that uses a semi-implicit 
Backward Euler to solve this system of equations and find values for a and b 
such that x and y oscillate in time.

For time-stepping using the semi-implicit method, evaluate the linear term 
in x at time t + ∆t, but evaluate the term x^2*y at time t. 
Hint for finding the oscillations: The initial concentrations of x and y 
need to be close to but not equal to the equilibrium concentrations. 
Determine the equilibrium values and use values that are close to this for 
the initial conditions.

Make a plot of x and y as functions of time showing a situation where you 
get oscillations, and one where the concentrations go to their equilibrium 
values by the end of the simulation. Make sure to tell what values you 
used for the initial conditions and for a and b for each case.
"""

# equilibrium when dx/dt = dy/dt = 0 => x = a + b, y = b/(a+b)^2 
def equilibrium(a, b):
    xeq = a + b
    yeq = b/((a + b)**2)
    return xeq, yeq

import numpy as np
import matplotlib.pyplot as plt

# define 2 given functions
def funcX(x, y, a):
    return  a-x + x**2*y

def funcY(x, y, b):
    return b - x**2*y

# apply semi-implicit Backward Euler method
def BackwardEuler (x0 , y0 , a, b, T, dt):
    time = np.arange(0, T, dt)
    xvals = np.zeros(len(time))
    xvals[0]  = x0
    yvals = np.zeros(len(time))
    yvals[0] = y0
    
    for i in range (1, len(time)):
        # evaluate a-x with time t+dt and x^2*y with time t
        xvals[i] = (xvals[i-1] + dt* (a+ xvals[i-1]**2 * yvals[i-1]))/(1 + dt)
        yvals[i] = yvals[i-1] + dt *funcY(xvals[i-1], yvals[i-1], b)
    
    return time, xvals, yvals

# define initial cons and solve for stable and unstable cases
T = 100
dt = 0.01

# larger b => unstable/oscillation. smaller b => stable/equillibrium
# yeq >> xeq => overshooting ==> unstable (like problem 1)

a_unstable = 0.01 
b_unstable = 1 

a_stable = 1.05   
b_stable = 1.05   

xeq_osc, yeq_osc = equilibrium(a_unstable, b_unstable)
xeq_stab, yeq_stab = equilibrium(a_stable, b_stable)

# using the hint: start close but not equal to equilibrium

x0_osc = xeq_osc + 0.1
y0_osc = yeq_osc + 0.1

x0_stab = xeq_stab + 0.05
y0_stab = yeq_stab + 0.05

time_osc, x_osc, y_osc = BackwardEuler(x0_osc, y0_osc, a_unstable, b_unstable, T, dt)
time_stab, x_stab, y_stab = BackwardEuler(x0_stab, y0_stab, a_stable, b_stable, T, dt)

# plot time

plt.subplot(1, 2, 1)
plt.plot(time_osc, x_osc, label='x (unstable)')
plt.plot(time_osc, y_osc, label='y (unstable)')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Oscillattion')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(time_stab, x_stab, label='x (stable)')
plt.plot(time_stab, y_stab, label='y (stable)')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Equilibrium')
plt.legend()

plt.tight_layout()
plt.show()




