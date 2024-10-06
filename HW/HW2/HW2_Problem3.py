# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 20:36:36 2024

@author: kietb
"""

"""
3. Write code that uses 4th order Runge-Kutta to solve
    dm/dt = c(1-m) - 10wm
    dc/dt = 5m(1-c) - 1.25c + S(t)
    dw/dt = 0.1(1-w) - 4mw

where S(t) = α when t is between 3 and 3.2, and is zero otherwise. 
Use initial conditions, m(0) = 0.0114, c(0) = 0.0090, and w(0) = 0.9374. 
Vary α to find a value for which the equilibrium switches to a new location. 
You should use a ∆t = 0.002 and run for a total time of 30. Make a
plot showing m, c, and w as functions of time when the equilibrium position 
switches. 
"""

import numpy as np
import matplotlib.pyplot as plt

# Define the system of equations
def equations(t, m, c, w, alpha):
    # conditions for S(t)
    if 3 <= t <= 3.2:
        S = alpha
    else:
        S = 0

    dm_dt = c * (1 - m) - 10 * w * m
    dc_dt = 5 * m * (1 - c) - 1.25 * c + S
    dw_dt = 0.1 * (1 - w) - 4 * m * w

    return np.array([dm_dt, dc_dt, dw_dt])

# apply the Fourth-order Runge-Kutta method
def RK4(f, y0, dt, T , alpha):
    time = np.arange(0, T, dt)
    mvals = np.zeros(len(time))
    cvals = np.zeros(len(time))
    wvals = np.zeros(len(time))

    
    mvals[0], cvals[0], wvals[0] = y0

    
    for i in range(1, len(time)):
        t = time[i - 1]
        y = np.array([mvals[i - 1], cvals[i - 1], wvals[i - 1]])

        k1 = dt * f(t, *y, alpha)
        k2 = dt * f(t + dt / 2, *(y + k1 / 2), alpha)
        k3 = dt * f(t + dt / 2, *(y + k2 / 2), alpha)
        k4 = dt * f(t + dt, *(y + k3), alpha)

        y_final = y + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)

        mvals[i], cvals[i], wvals[i] = y_final

    return time, mvals, cvals, wvals

# define initial cons
m0 = 0.0114
c0 = 0.0090
w0 = 0.9374
y0 = [m0, c0, w0]

T = 30
dt = 0.002
alphas = np.linspace(0, 20, 100)

for alpha in alphas:
   time, mvals, cvals, wvals = RK4(equations, y0, dt, T, alpha)

# plot
plt.plot(time, mvals, label="m(t)")
plt.plot(time, cvals, label="c(t)")
plt.plot(time, wvals, label="w(t)")
plt.xlabel("Time")
plt.ylabel("Values")
plt.title("RK4 Method")
plt.legend()
plt.grid(True)
plt.show()