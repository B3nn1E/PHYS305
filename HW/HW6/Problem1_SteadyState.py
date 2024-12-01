# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 18:22:27 2024

@author: kietb
"""

'''
Problem 1: Determine the steady state profile of the following 
equation:
 ∂C / ∂t =  ∂^2C/ ∂x^2 - ∂/ ∂x(sin(x)C)
on a domain that goes from x = 0 to x = 4π using periodic boundary 
conditions. Make sure to use the appropriate upwind derivatives. 
'''
# Step 1: Divide N domain into grid points
# Step 2: Create Laplacian matrix for the second deriv
# Step 3: Create advection matrix fir the first deriv with sinx
# Step 4: Combine those 2 matrices
# Step 5: Normalization and solve!

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
def SteadyState(N = 100, L = 4 * np.pi):
    # Solve the steady state equation on the domain [0, 4pi]
    
    # Step 1: Define grid
    x, dx = np.linspace(0, L, N, endpoint = False, retstep = True)
    # endpoint to make sure 4pi does not include for bound cond
    sin_x = np.sin(x)
    
    # Step 2: Set up Laplacian matrix using central different
    main_diag = -2 * np.ones(N) # coeff of Ci for central diff is -2
    off_diag = np.ones(N - 1) # set C_[i+1] and C_[i-1] = 1
    lap = diags([off_diag, main_diag, off_diag], [-1, 0, 1], shape = (N, N))
    lap = lap / dx**2
    
    # Periodic boundary conditions
    # "Link" the first row (x = 0) to the last column(x = L)
    # "Link" the last row (x = L) to the first column (x = 0)
    lap = lap.tolil() # make matrix to list of list for easy modification
    lap[0, -1] = 1 / dx**2 # top right to 1 / dx**2
    lap[-1, 0] = 1 / dx**2 # bottom left to 1 / dx**2
    
    # Step 3: Set up advection matrix 
    adv = diags([-sin_x, sin_x], [-1, 0], shape = (N, N))
    adv = adv / dx
    
    # Periodic boundary conditions
    # Similar process to the Laplacian
    adv = adv.tolil()
    adv[0, -1] = -sin_x[-1] / dx # take sin(x) at the last grid point
    adv[-1, 0] = sin_x[0] / dx # take sin(x) at the first grid point
    
    # Step 4: Combine Laplacian and advection for matrix A
    A = lap - adv
    
    # Step 5: Now solve for Ax = 0 (steady-state)
    # AC = b so b = 0
    # Normalization first by setting the first element of A = 1
    # Needed to avoide null space, ensuring unique solution
    b = np.zeros(N)
    b[0] = 1  # Normalization condition
    A[0, 0] = 1  # Ensure unique solution
    
    C = spsolve(A, b)
    
    # Plot the steady-state solution
    plt.plot(x, C, label = 'Steady-state')
    plt.xlabel('x')
    plt.ylabel('C(x)')
    plt.title('Steady-state Solution')
    plt.grid(True)
    plt.legend()
    plt.show()
    
    return x, C

# Run the function
x, C = SteadyState()
    
    