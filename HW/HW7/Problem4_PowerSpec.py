# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:43:32 2024

@author: kietb
"""

'''
Problem 4:
    Download the data file Problem4_data.txt. The first column of data 
    in this file is time, and the second column is the amplitude of a wave. 
    Determine the power spectrum of the data and tell what the three most 
    dominant frequencies are. 
'''

import numpy as np
import matplotlib.pyplot as plt

def PowerSpectrum(datafile):

    
    # Load data from file
    file = open(datafile, 'r')
    data = []
    for line in file:
        data.append(line.split())
    file.close()

    # Convert data to numpy array
    data = np.array(data, dtype=float)
    time = data[:, 0]
    amplitude = data[:, 1]

    # Perform FFT
    fft_result = np.fft.fft(amplitude)
    N = len(amplitude)
    freq = np.fft.fftfreq(N, d=(time[1] - time[0]))

    # Compute power spectrum
    power_spectrum = np.abs(fft_result)**2

    # Plot the power spectrum
    plt.figure(figsize=(8, 6))
    plt.plot(freq[:N // 2], power_spectrum[:N // 2], label='Power Spectrum')
    plt.title('Power Spectrum')
    plt.xlabel('Frequency')
    plt.ylabel('Power')
    plt.grid()
    plt.legend()
    plt.show()

    return dominant_frequencies

# Load file and apply function
dominant_frequencies = PowerSpectrum(r'C:\Users\kietb\OneDrive\Desktop\Suffering\Undergrad\PHYS305\HW\HW7\Problem4_data.txt')
print("Three most dominant frequencies:", dominant_frequencies)
