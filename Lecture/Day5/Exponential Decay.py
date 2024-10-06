# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 15:32:45 2024

@author: kietb
"""

# N0 = 1, k = 1, dt = 0.1

import numpy as np
from matplotlib.pyplot import plot

time, dt = np.linspace(0, 10, 100, 'retstep', 'true', 'float')

Nf = np.zeros((100,1), 'float') # forward method
Nb = np.zeros((100,1), 'float') # backward method

# initial con
Nf[0] = 1
Nb[0] = 1

for i in range (1, 20):
    Nf[i] = (1-dt) * Nf[i-1]
    Nb[i] = Nb[i-1]/(1+dt)
    
plot (time, Nf, time, Nb)

# dt = 0.5
time, dt = np.linspace(0, 10, 20, 'retstep', 'true', 'float')

Nf = np.zeros((20,1), 'float') # forward method
Nb = np.zeros((20,1), 'float') # backward method

# initial con
Nf[0] = 1
Nb[0] = 1

for i in range (1, 20):
    Nf[i] = (1-dt) * Nf[i-1]
    Nb[i] = Nb[i-1]/(1+dt)
    
plot (time, Nf, time, Nb)

# dt = 1
time, dt = np.linspace(0, 10, 10, 'retstep', 'true', 'float')

Nf = np.zeros((10,1), 'float') # forward method
Nb = np.zeros((10,1), 'float') # backward method

# initial con
Nf[0] = 1
Nb[0] = 1

for i in range (1, 10):
    Nf[i] = (1-dt) * Nf[i-1]
    Nb[i] = Nb[i-1]/(1+dt)
    
plot (time, Nf, time, Nb)

# dt = 2
time, dt = np.linspace(0, 10, 5, 'retstep', 'true', 'float')

Nf = np.zeros((5,1), 'float') # forward method
Nb = np.zeros((5,1), 'float') # backward method

# initial con
Nf[0] = 1
Nb[0] = 1

for i in range (1, 5):
    Nf[i] = (1-dt) * Nf[i-1]
    Nb[i] = Nb[i-1]/(1+dt)
    
plot (time, Nf, time, Nb)
