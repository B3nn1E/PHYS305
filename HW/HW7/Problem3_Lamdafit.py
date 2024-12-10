# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:01:42 2024

@author: kietb
"""
'''
Problem 3:
    Using the same data file, fit the data using
            λ * ∂^2f / ∂x^2 - f = f_data
    How does your fit depend on λ? Is there a choice for λ that 
    you feel works well? Explain why.
    
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded 
# use to solve banded matrix, quicker to use this function


def LambdaFit(datafile, lambdas):
    
    # Load data
    file = open(datafile, 'r')
    data = []
    for line in file:
        data.append(line.split())
    file.close()

    data = np.array(data, dtype='float')
    x = data[:, 0]
    y = data[:, 1]
    n = len(x)
    dx = x[1] - x[0]

    fits = {}

    plt.figure(figsize=(12, 8))
    # plt.scatter(x, y, label='Data', color='blue', s=10)
    plt.plot(x, y, label='Data', color='blue')

    for lambd in lambdas:
        # Define matrix
        main_diag = (2 + 2 * lambd / dx**2) * np.ones(n)
        off_diag = -lambd / dx**2 * np.ones(n - 1)

        # Banded matrix 
        # Used usually when discretizing differential equations 
        ab = np.zeros((3, n))
        ab[0, 1:] = off_diag  # Upper diag
        ab[1, :] = main_diag  # Main diag
        ab[2, :-1] = off_diag  # Lower diag

        # Solve for smoothed fit
        # Large 𝜆 fit is smoother but less accurate
        b = y.copy()
        smooth_y = solve_banded((1, 1), ab, b)
        fits[lambd] = smooth_y

        # Plot fit
        plt.plot(x, smooth_y, label=f'Lambda = {lambd}')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Lambda Fit')
    plt.legend()
    plt.grid()
    plt.show()

    return fits

# Define lambdas' values
lambdas = [0.001, 0.01, 0.1, 1, 10]
fits = LambdaFit(r'C:\Users\kietb\OneDrive\Desktop\Suffering\Undergrad\PHYS305\HW\HW7\Problem2_data.txt', lambdas)
