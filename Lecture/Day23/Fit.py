# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 16:47:29 2024

@author: kietb
"""

import numpy as np
import matplotlib.pyplot as plt

def FitDataLinear(datafile):
    file = open(datafile, 'r')
    data = []
    for line in file:
        data.append(line.split())

    data = np.array(data, dtype='float')
    N, d = data.shape

    a1 = (N * np.sum(data[:, 1] * data[:, 0]) - np.sum(data[:, 1]) * np.sum(data[:, 0])) / \
         (N * np.sum(data[:, 0]**2) - np.sum(data[:, 0])**2)
    a0 = (np.sum(data[:, 1]) - a1 * np.sum(data[:, 0])) / N

    x = np.linspace(data[0, 0], data[-1, 0], 100)
    y = a0 + a1 * x

    plt.plot(data[:, 0], data[:, 1], 'b')
    plt.plot(x, y, 'r')
    plt.title('Linear Fit')
    plt.show()

    return a0, a1

def FitDataCubic(datafile):
    file = open(datafile, 'r')
    data = []
    for line in file:
        data.append(line.split())

    data = np.array(data, dtype='float')
    N, d = data.shape
    x = data[:, 0]
    y = data[:, 1]

    M = np.zeros((4, 4))
    M[:, 0] = [N, np.sum(x), np.sum(x**2), np.sum(x**3)]
    M[:, 1] = [np.sum(x), np.sum(x**2), np.sum(x**3), np.sum(x**4)]
    M[:, 2] = [np.sum(x**2), np.sum(x**3), np.sum(x**4), np.sum(x**5)]
    M[:, 3] = [np.sum(x**3), np.sum(x**4), np.sum(x**5), np.sum(x**6)]

    b = [np.sum(y), np.sum(x * y), np.sum(x**2 * y), np.sum(x**3 * y)]
    a = np.linalg.solve(M, b)

    y_fit = a[0] + a[1] * x + a[2] * x**2 + a[3] * x**3

    plt.plot(data[:, 0], data[:, 1], 'b')
    plt.plot(x, y_fit, 'r')
    plt.title('Cubic Fit')
    plt.show()

    return a

def ComputePowerSpectrum(datafile):
    file = open(datafile, 'r')
    data = []
    for line in file:
        data.append(line.split())

    data = np.array(data, dtype=float)
    time = data[:, 0]
    amplitude = data[:, 1]

    # FFT
    f_fft = np.fft.fft(amplitude)
    N = len(amplitude)
    freq = np.fft.fftfreq(N, d=(time[1] - time[0]))

    # Compute power spectrum
    power = np.abs(f_fft)**2

    # Plot power spectrum
    plt.plot(freq[:N // 2], power[:N // 2])
    plt.title('Power Spectrum')
    plt.xlabel('Frequency')
    plt.ylabel('Power')
    plt.grid()
    plt.show()
