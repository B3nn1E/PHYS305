# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:33:56 2024

@author: kietb
"""

"""
import math
v = [1+i*0.2 for i in range(30)]
logV = math.log(v) #math can't handle list
"""
import numpy as np
v = np.array([1+i*0.2 for i in range(30)], 'float')
logV = np.log(v)
logV

v, dv = np.linspace(1,5.9,30,'retstep', 'true', 'float')
v
dv

