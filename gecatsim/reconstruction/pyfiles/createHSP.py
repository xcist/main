# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import numpy as np
import math
from scipy.interpolate import interp1d

def createHSP(Length, kernelType):
    HS = np.zeros(Length)
    Center = int((Length) / 2)
    PI = 3.14159265358979
    nn = Length
    nn2 = nn * 2

    if kernelType == "R-L":
        HS[0] = 0
        HS[Center] = 0.25
        for i in range(1, Center):
            HS[i] = -np.power((math.sin(PI * (i - Center) / 2)), 2) / (PI * PI * (i - Center) * (i - Center))

        for k in range((Center + 1), Length):
            HS[k] = -np.power((math.sin(PI * (k - Center) / 2)), 2) / (PI * PI * (k - Center) * (k - Center))
        k = int(nn / 2)
        TempF = np.zeros(nn2)
        TempF[0:k] = HS[k:nn]
        TempF[k + nn:nn2] = HS[0:k]
        HS = TempF * complex(0, 1)
        FFT_F = np.fft.fft(HS)

    elif kernelType == "S-L":
        for i in range(Length):
            HS[i] = -2 / (PI * PI * (4 * (i - Center) * (i - Center) - 1))
        k = int(nn / 2)
        TempF = np.zeros(nn2)
        TempF[0:k] = HS[k:nn]
        TempF[k + nn:nn2] = HS[0:k]
        HS = TempF * complex(0, 1)
        FFT_F = np.fft.fft(HS)

    elif kernelType == "soft":
        x = np.array([0, 0.25, 0.5, 0.75, 1])
        y = np.array([1, 0.815, 0.4564, 0.1636, 0])
        # y = np.array([1, 1.0485, 1.17, 1.2202, 0.9201])
        # y= np.array([1, 0.9338, 0.7441, 0.4425, 0.0531])
        f = interp1d(x, y, kind='quadratic')
        FFT_F = np.zeros(nn2)
        for i in range(nn):
            FFT_F[i] = f((i)/nn)*0.997 * (i+0.003) / nn2
            FFT_F[nn2 - i - 1] = f((i)/nn)* 0.997*(i + 1 + 0.003) / nn2
        FFT_F= FFT_F * complex(0,1)

    elif kernelType == "standard":
        x = np.array([0, 0.25, 0.5, 0.75, 1])
        # y = np.array([1, 0.815, 0.4564, 0.1636, 0])
        # y = np.array([1, 1.0485, 1.17, 1.2202, 0.9201])
        y = np.array([1, 0.9338, 0.7441, 0.4425, 0.0531])
        f = interp1d(x, y, kind='quadratic')
        FFT_F = np.zeros(nn2)
        for i in range(nn):
            FFT_F[i] = f((i) / nn) * 0.997 * (i + 0.003) / nn2
            FFT_F[nn2 - i - 1] = f((i) / nn) * 0.997 * (i + 1 + 0.003) / nn2
        FFT_F = FFT_F * complex(0, 1)

    elif kernelType == "bone":
        x = np.array([0, 0.25, 0.5, 0.75, 1])
        # y = np.array([1, 0.815, 0.4564, 0.1636, 0])
        y = np.array([1, 1.0485, 1.17, 1.2202, 0.9201])
        # y = np.array([1, 0.9338, 0.7441, 0.4425, 0.0531])
        f = interp1d(x, y, kind='quadratic')
        FFT_F = np.zeros(nn2)
        for i in range(nn):
            FFT_F[i] = f((i) / nn) * 0.997 * (i + 0.003) / nn2
            FFT_F[nn2 - i - 1] = f((i) / nn) * 0.997 * (i + 1 + 0.003) / nn2
        FFT_F = FFT_F * complex(0, 1)

    else: 
        raise Exception("******** Error! An unsupported kernel was specified: {:s}. ********".format(kernelType))

    return FFT_F