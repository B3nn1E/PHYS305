# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 15:32:08 2024

@author: kietb
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def double_deriv(dx, C):
    D = 0.5
    dCd2 = np.zeros(len(C))
    for i in range(1, len(C) - 1):
        dCd2[i] = D * (C[i + 1] - 2 * C[i] + C[i - 1]) / (dx * dx)
    return dCd2

def deriv(dx, C):
    v = 0.25
    dC = np.zeros(len(C))
    dC[0] = -1 + 0.5 * C[0]
    dC[1:-1] = v * (C[2:] - C[:-2]) / (2 * dx)
    dC[-1] = 1 + 0.5 * C[-1]
    return dC

def main(dx, L, dt, t_f):
    N_steps = int(t_f / dt)
    N = int(L / dx)
    C = np.ones((N, N_steps))
    for i in range(1, N_steps):
        dC = deriv(dx, C[:, i - 1])
        dC2 = double_deriv(dx, C[:, i - 1])
        C[:, i] = C[:, i - 1] + dt * (0.5 * dC - 0.25 * dC2 - C[:, i - 1])
    return C

L = 100
dx = 1
t_f = 10
dt = 0.1

C = main(dx, L, dt, t_f)

fig, ax = plt.subplots()
line, = ax.plot(np.linspace(0, L, int(L / dx)), C[:, 0])

def update(frame):
    line.set_ydata(C[:, frame])
    ax.set_title(f'Time: {frame * dt:.2f} s')
    return line,

ani = FuncAnimation(fig, update, frames=range(0, C.shape[1], 5), blit=True)
plt.show()
