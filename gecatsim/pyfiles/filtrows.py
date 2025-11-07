# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

# -----------------------------------------------------------------------
# Aim
#   Function to convolve the rows of x with the 1-D filter h
#   The filter must already be the correct length, and in FFT order
# -----------------------------------------------------------------------

import numpy as np

def filtrows(x, h):
    X = np.fft.fft(x, len(h), axis=1)
    H = np.ones((x.shape[0], 1)) * np.fft.fft(h)
    Y = X * H
    y = np.real(np.fft.ifft(Y, len(h), axis=1))
    y = y[:, :x.shape[1]]
    return y