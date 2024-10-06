# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 17:36:21 2024

@author: kietb
"""

"""
Problem 4:
Write code that integrates sin x between two user-defined limits
using the lefthand, righthand and trapezoidal rules. 
Determine the error in these integrations as a function of the 
step size
"""
import math
def sinLH (a,b,Num):
    
# define step size
    dF = (b-a)/(Num-1)

    first = a # start on the left side/first sub interval
    sin_left = 0
    
    for i in range(Num-1):
        
        sin_left += math.sin(first)
        first += dF

    Error_left = abs(sin_left - (math.cos(a)-math.cos(b)))
    return (sin_left, Error_left)


def sinRH (a,b,Num):
    
# define step size
    dF = (b-a)/(Num-1)

    first = a + dF # start on the left side/first sub interval
    sin_right = 0
    
    for i in range(Num-1):
        
        sin_right += math.sin(first)
        first += dF

    Error_right = abs(sin_right - (math.cos(a)-math.cos(b)))
    return (sin_right, Error_right)



def sinTrap (a,b,Num):

# define step size
    dF = (b-a)/(Num-1)  
    
# initialize Volume and Work for the trapezoidal rule
    sin_trap = 0
    first = a # left hand first sub interval
    
# calculate Work done (Trapezoidal)
    for i in range(Num-1):
        second = first + dF #right hand second sub interval
        # Trapezoidal is the average of left and right
        sin_trap += 0.5 * (math.sin(first) + math.sin(second))
        first = second
        
    Error_trap = abs(sin_trap - (math.cos(a)-math.cos(b)))
    return (sin_trap, Error_trap)


# set Num in increments of 10 
Num = [10 + 10*i for i in range (10)]
Sin_Left= [0 for i in range(10)]
Error_Left = [0 for i in range(10)]
Sin_Right= [0 for i in range(10)]
Error_Right = [0 for i in range(10)]
Sin_Trap= [0 for i in range(10)]
Error_Trap = [0 for i in range(10)]

for i in range(10):
    Sin_Left[i], Error_Left[i] = sinLH(0, 3.14, Num[i])
    Sin_Right[i], Error_Right[i] = sinRH(0,3.14,Num[i])
    Sin_Trap[i], Error_Trap[i] = sinTrap(0,3.14,Num[i])

import matplotlib.pyplot as plt
plt.plot(Num, Error_Left, 'o', label='LH')
plt.plot(Num, Error_Right, 'o', label='RH')
plt.plot(Num, Error_Trap, 'o', label='Trap')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Num')
plt.ylabel('Error')
plt.legend()
plt.show()