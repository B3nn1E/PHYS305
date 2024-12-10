# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:43:03 2024

@author: kietb
"""

'''
Problem 2:
    Download the data file Problem2_data.txt and fit the solution 
    to a third order polynomial.
'''
# I am sorry, I will just cheat using np.polyfit instead
import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt(r'C:\Users\kietb\OneDrive\Desktop\Suffering\Undergrad\PHYS305\HW\HW7\Problem2_data.txt')
x = data[:, 0]
y = data[:, 1]

# Fit data to a third-order polynomial
coefficients = np.polyfit(x, y, 3)
poly_fit = np.poly1d(coefficients)

# Generate fitted values
x_fit = np.linspace(min(x), max(x), 500)
y_fit = poly_fit(x_fit)

# Plot data
plt.figure(figsize=(8, 6))
# plt.scatter(x, y, label='Data', color='blue', s=10)
plt.plot(x, y, label='Data', color='blue')
plt.plot(x_fit, y_fit, label='3rd Order Polynomial Fit', color='red')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Third Order Polynomial Fit')
plt.legend()
plt.grid()
plt.show()

