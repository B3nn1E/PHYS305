# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 15:32:11 2024

@author: kietb
"""

from Numerical_integration import IdealGasWorkLH

# set Num increments of 10
Num = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# better way
Num = [10 + 10*i for i in range (10)]
Work = [0 for i in range(10)]
Error = [0 for i in range(10)]

for i in range(10):
    Work[i], Error[i] = IdealGasWorkLH(1,10,Num[i])
    
    
from matplotlib.pyplot import plot 
plot (Num, Error) # line plot
plot (Num, Error, 'o') # scatter plot
plot (Num, Error, 'or') #red scatter plot

import math 
plot(math.log(Num), math.log(Error), 'o') # can't bc Num and Error are lists

logN = [0 for i in range (10)]
logE = [0 for i in range(10)]

for i in range(10):
    logN[i] = math.log(Num[i])
    logE[i] = math.log(Error[i])
    
plot(logN, logE, 'or')