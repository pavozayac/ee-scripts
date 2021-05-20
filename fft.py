import math
import sys
import numpy as np

def fft (*coefficients) -> list:
    points = []

    n = len(coefficients)

    if n == 1:
        return [coefficients[0]]
    else:
        even = fft(*coefficients[0::2])
        odd = fft(*coefficients[1::2])

        print(n/2)
        print(int(n/2))
        for k in range(int(n/2)):
            #print(k)
            runity: complex = math.e**(-2*math.pi*1j*k/n)

            points.insert(k, even[k] + runity*odd[k])
            points.insert(int(k+n/2), even[k] - runity*odd[k])

    return points



def DFT_slow(x):
    """Compute the discrete Fourier Transform of the 1D array x"""
    x = np.asarray(x, dtype=float)
    N = x.shape[0]
    n = np.arange(N)
    k = n.reshape((N, 1))
    M = np.exp(-2j * np.pi * k * n / N)
    return np.dot(M, x)

def FFT(x):
    """A recursive implementation of the 1D Cooley-Tukey FFT"""
    x = np.asarray(x, dtype=float)
    N = x.shape[0]
    
    if N % 2 > 0:
        raise ValueError("size of x must be a power of 2")
    elif N <= 32:  # this cutoff should be optimized
        return DFT_slow(x)
    else:
        X_even = FFT(x[::2])
        X_odd = FFT(x[1::2])
        factor = np.exp(2j * np.pi * np.arange(N) / N)
        return np.concatenate([X_even + factor[:N / 2] * X_odd,
                               X_even + factor[N / 2:] * X_odd])

if __name__ == '__main__':
    sys.setrecursionlimit(2000)

    print(np.array(fft(1, 2, 3, 4))*np.array(fft(2, 1, 4, 3)))