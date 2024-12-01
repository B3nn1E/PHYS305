# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 16:16:40 2024

@author: kietb
"""

import numpy as np
from time import time
from Gauss_Jordan import GJ
from scipy.sparse import csr_matrix
from scipy.linalg import spsolve
from matplotlib.pyplot import plot 

def GJ2(A,b):
    # performs Gauss-Jordan elimination to solve A x = b 
    
    # derermine the shape of A
    W, L = A.shape
    
    # define matrix B so that A does not get overwritten
    B = np.array(A, dtype = 'float')
    
    # initialize array for solution x
    x = np.array(b, dtype = 'float')
    
    vals = np.array(k for k in range(W))

    # begin iteration
    for j in range(W):
        
    # check if diagonal element is zero, pivot if it is
        if B[j, j] == 0:
            
            # in the jth column find number points that are nonzero
            nonzero = [k for k, e in enumerate(B[j:,j]) if e != 0]
            
            # grab the number of the first nonzero element
            val = nonzero[0]
            
            # exact the current row and first row with nonzero diagonal
            b1 = B[j + val,:].copy()
            b2 = B[j,:].copy()
            
            # swap those rows
            B[j,:] = b1
            B[j + val,:] = b2
            
            # get the current row and first row with nonzero diagonal
            c1 = x[j + val]
            c2 = x[j]
            
            # swap those rows
            x[j] = c1
            x[j + val] = c2
            
        # divide the jth row of  B by  the diagonal element, likewise for q
        norm = B[j, j] 
        B[j,:] = B[j,:] / norm
        x[j] = x[j] / norm
        
        # for all i not equal to j subtract correct multiple of jth row to
        # # get zero in the jth column of the ith row
        
        for i in vals[(vals != j) & (B[:,j] != 0)]:
            norm = B[i,j]
            B[i,:] = B[i,:] - norm * B[j,:]
            x[i] = x[i] - norm * x[j]
            
            
        return x

N = 500

# define a large second deriv matrix
A = np.zeros((N,N))
A[0,0] = 1
for i in range(1, N - 1):
    A[i, i - 1:i + 2] = [1, -2, -1]
A [N - 1, N - 1] = 1

# define the solution vector
b = np.ones(N)

# put in correct boundary conditions
b[0] = 0
b[N - 1] = 0

# for GJ   
t1 = time()
x1 = GJ(A, b)
t2 = time()
print('It took{:} s to solve using GJ.'.format(t2 - t1))

# for GJ2
t3 = time()
x2 = GJ2(A, b)
t4 = time()
print('It took{:} s to solve using GJ2.'.format(t4 - t3))

# for GJ2
t5 = time()
x3 = np.alg.solve(A, b)
t6 = time()
print('It took{:} s to solve using Numpy.'.format(t5 - t6))

# for GJ2
A_sparse = csr_matrix(A)
t7 = time()
x4 = spsolve(A_sparse, b)
t8 = time()
print('It took{:} s to solve using sparse matrices'.format(t7 - t8))

plot(x1, 'b')
plot(x2, 'r')
plot(x3, 'g')
plot(x4, 'k')