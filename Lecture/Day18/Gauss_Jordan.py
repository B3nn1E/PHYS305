# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 15:54:12 2024

@author: kietb
"""

import numpy as np

def GJ(A,b):
    # performs Gauss-Jordan elimination to solve A x = b 
    
    # derermine the shape of A
    W, L = A.shape
    
    # define matrix B so that A does not get overwritten
    B = np.array(A, dtype = 'float')
    
    # initialize array for solution x
    x = np.array(b, dtype = 'float')
    
    # begin iteration
    for j in range(W):
        
    # check if diagonal element is zero, pivot if it is
        if B[j, j] == 0 :
    
            # in the jth column find number points that are nonzero
            nonzero = [k for k, e in enumerate(B[j:,j]) if e != 0]
            
            # grab the number of the first nonzero element
            val = nonzero[0]
            
            # exact the current row and first row with nonzero diagonal
            b1 = np.copy(B[j - val,:])
            b2 = np.copy(B[j,:])
            
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
        
        # for all i not equal to j
        vals = [k for k in range(W)]
        vals.remove(j)
        
        # subtract correct multiple of jth row to get zero in the jth
        # column of the ith row
        for i in vals:
            norm = B[i,j]
            vals.remove(j)
            
        # subtract correct multiple of jth row to get zero in the jth
        # column of the ith row
        for i in vals:
            norm = B[i,j]
            B[i,:] = B[i,:] - norm * B[j,:]
            x[i] = x[i] - norm * x[j]
            
        return x
    