# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:53:01 2024

@author: kietb
"""
import math 

def IdealGasWorkLH (Vi,Vf,Num):
    """    
    Computes the work done by an ideal gas wehn it expands from Vi to Vf 
    using the left hand rule
    """
# define parameters 
    N = 2.2e23 #For a typical canister holding 16g of CO2
    kbT = 4e-14 #ergs

# define step size
    dV = (Vf-Vi)/(Num-1)

# initialize Volume and Work
    V = Vi
    Work = 0
    
    """
    for and range command
    for i in range (2,10)
    will execute from i = 2 to i = 9
    """
    
    # calculate Work done
    for i in range(Num-1):
        
        Work = Work + (N*kbT/V)*dV
        V = V + dV
    
    Error = abs(Work - (N*kbT)*math.log(Vf/Vi))
    
    print ('The Work done by the gas when it expands from {:2.2f} ml to {:2.2f} ml is {:2.2f} ergs'.format(Vi,Vf,Work))
    return Work, V, Error
