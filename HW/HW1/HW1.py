# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 14:34:24 2024

HW1: Ben Phan
"""

"""
Problem 1:
Write code that takes a temperature in Fahrenheit that is 
input by the user from the keyboard and converts it to Celsius. 
The output should be print to the screen.
"""
f = float(input('Please input your temperature in Fahrenheit '))
c = print((f-32)*(5/9))

"""
Problem 2:
Write a function that computes the work done by an ideal gas when it
expands from one volume to another. The volumes and the number of 
steps should be arguments to the function. Compute the work done using
the lefthand, righthand, and trapezoidal rules. The code should 
return the three solutions
"""

import math 

def IdealGasWork (Vi,Vf,Num):
    """    
    Computes the work done by an ideal gas when it expands from Vi 
    to Vf using the left hand rule, right hand rule, and trapezoidal
    rule
    
    """
# define parameters 
    N = 2.2e23 #For a typical canister holding 16g of CO2
    kbT = 4e-14 #ergs

# define step size
    dV = (Vf-Vi)/(Num-1)

# initialize Volume and Work for the left hand rule
    V_left = Vi # start on the left side/first sub interval
    Work_left = 0
    
# calculate Work done (LH)
    for i in range(Num-1):
        
        Work_left += (N*kbT/V_left)*dV
        V_left += + dV
        
# initialize Volume and Work for the right hand rule
    V_right = Vi + dV # start on the right side/second sub interval
    Work_right = 0
    
# calculate Work done (RH)
    for i in range(Num-1):
        
        Work_right += (N*kbT/V_right)*dV
        V_right += dV
        
# initialize Volume and Work for the trapezoidal rule
    Work_trap = 0
    V_first = Vi # left hand first sub interval
    
# calculate Work done (Trapezoidal)
    for i in range(Num-1):
        V_second = V_first + dV #right hand second sub interval
        # Trapezoidal is the average of left and right
        Work_trap += 0.5 * (N*kbT/V_first + N*kbT/V_second)*dV
        V_first = V_second
        
    return (Work_right, Work_left, Work_trap)

"""
Problem 3: 
Write code that uses the function from Problem 3 to compute the 
error between each of the three integration methods and the analytic
solution. Use at least 7 different step sizes for each
method and determine how the error scales with step size.
"""

def IdealGasWorkLH (Vi,Vf,Num): # break down function 2 to LH, RH, Trap 
    
# define parameters 
    N = 2.2e23 #For a typical canister holding 16g of CO2
    kbT = 4e-14 #ergs

# define step size
    dV = (Vf-Vi)/(Num-1)

# initialize Volume and Work for the left hand rule
    V_left = Vi # start on the left side/first sub interval
    Work_left = 0
    
# calculate Work done (LH)
    for i in range(Num-1):
        
        Work_left += (N*kbT/V_left)*dV
        V_left += dV
        
# calculate error (LH)
    Error_left = abs(Work_left - (N*kbT)*math.log(Vf/Vi))
    return (Work_left, Error_left)


def IdealGasWorkRH (Vi,Vf,Num):
    
# define parameters 
    N = 2.2e23 #For a typical canister holding 16g of CO2
    kbT = 4e-14 #ergs

# define step size
    dV = (Vf-Vi)/(Num-1)   
     
# initialize Volume and Work for the right hand rule
    V_right = Vi + dV # start on the right side/second sub interval
    Work_right = 0
    
# calculate Work done (RH)
    for i in range(Num-1):
        
        Work_right += (N*kbT/V_right)*dV
        V_right += dV
    
# calculate error (RH)
    Error_right = abs(Work_right - (N*kbT)*math.log(Vf/Vi))
    return (Work_right, Error_right)


def IdealGasWorkTrap (Vi,Vf,Num):
    
# define parameters 
    N = 2.2e23 #For a typical canister holding 16g of CO2
    kbT = 4e-14 #ergs

# define step size
    dV = (Vf-Vi)/(Num-1)  
    
# initialize Volume and Work for the trapezoidal rule
    Work_trap = 0
    V_first = Vi # left hand first sub interval
    
# calculate Work done (Trapezoidal)
    for i in range(Num-1):
        V_second = V_first + dV #right hand second sub interval
        # Trapezoidal is the average of left and right
        Work_trap += 0.5 * (N*kbT/V_first + N*kbT/V_second)*dV
        V_first = V_second
        
# calculate error (Trap) 
    Error_trap = abs(Work_trap - (N*kbT)*math.log(Vf/Vi))
    return (Work_trap, Error_trap)


# set Num in increments of 10 
Num = [10 + 10*i for i in range (10)]
Work_Left= [0 for i in range(10)]
Error_Left = [0 for i in range(10)]
Work_Right= [0 for i in range(10)]
Error_Right = [0 for i in range(10)]
Work_Trap= [0 for i in range(10)]
Error_Trap = [0 for i in range(10)]

for i in range(10):
    Work_Left[i], Error_Left[i] = IdealGasWorkLH(1,10,Num[i])
    Work_Right[i], Error_Right[i] = IdealGasWorkRH(1,10,Num[i])
    Work_Trap[i], Error_Trap[i] = IdealGasWorkTrap(1,10,Num[i])

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

"""
Problem 4:
Write code that integrates sin x between two user-defined limits
using the lefthand, righthand and trapezoidal rules. 
Determine the error in these integrations as a function of the 
step size
"""
def sinLH (a,b,Num):
    
# define step size
    dF = (b-a)/(Num-1)
    
# initialize
    first = a # start on the left side/first sub interval
    sin_left = 0

# calculate sin (x) (LH)
    for i in range(Num-1):
        
        sin_left += math.sin(first)
        first += dF
 
# multiply by step size  
    sin_left *= dF
    Error_left = abs(sin_left - (math.cos(a)-math.cos(b)))
    return (sin_left, Error_left)


def sinRH (a,b,Num):
    
# define step size
    dF = (b-a)/(Num-1)

# initialize
    first = a + dF # start on the left side/first sub interval
    sin_right = 0
 
# calculate sin (x) (RH)
    for i in range(Num-1):
        
        sin_right += math.sin(first)
        first += dF

# multiply by step size  
    sin_right *= dF
    Error_right = abs(sin_right - (math.cos(a)-math.cos(b)))
    return (sin_right, Error_right)



def sinTrap (a,b,Num):

# define step size
    dF = (b-a)/(Num-1)  
    
# initialize 
    sin_trap = 0
    first = a # left hand first sub interval
    
# calculate sin(x) (Trapezoidal)
    for i in range(Num-1):
        second = first + dF #right hand second sub interval
        # Trapezoidal is the average of left and right
        sin_trap += 0.5 * (math.sin(first) + math.sin(second))
        first = second
        
# multiply by step size   
    sin_trap *= dF
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
plt.plot(Num, Error_Trap, 'o', label='Trapo')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Num')
plt.ylabel('Error')
plt.legend()
plt.show()

"""
Problem 5e: 
Numerically compute the integral and tell what value you find. 
Tell what values you used for λ and ∆x
"""
# write down the function 
def func(x):
    if x == 0: # to avoid divide by 0 since the  lower limit is 0 
        return 0  # function approach 0 as x = 0
    return (1 - math.exp(-x)) / (x * (1 + x**2))


def funcTrap (a,b,Num): # I choose trapezoidal approx for this 

# define step size
    dF = (b-a)/(Num-1)  
    
    func_trap = 0
    first = a # left hand first sub interval
    
# calculate the integral in trapezoidal 
    for i in range(Num-1):
        second = first + dF #right hand second sub interval
        # Trapezoidal is the average of left and right
        func_trap += 0.5 * (func(first) + func(second))
        first = second
        
 # multiply by step size     
    func_trap *= dF
    return (func_trap)

# set the upper, lower limits and steps
# integer since 1e is considered float
a = int(1e-8)  
b = 1000  
Num = int(1e8)

result = funcTrap(a, b, Num)
print(result)
