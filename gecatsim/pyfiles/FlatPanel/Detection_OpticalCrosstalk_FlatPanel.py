# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import numpy as np
import matplotlib.pyplot as plt

def Detection_OpticalCrosstalk_FlatPanel(viewIn, cfg):
    '''
    Models the optical crosstalk of the flat panel detector. 
    It applies a Lorentzian in frequency domain (which is characterized by a fitting parameter scanner.FlatPanel_H) 
    to the measured view after photon absorption and before electronic gain conversion

    Author: Junyuan Li, Department of Biomedical Engineering, Johns Hopkins University
    '''
    viewOut_size = int(cfg.scanner.detectorColCount*cfg.physics.FlatPanel_OSfactor)
    viewOut = viewIn.reshape(viewOut_size, viewOut_size)
    
    
#   Model the optical crosstalk as a Lorentzian function in frequency domain
    lorentz_row_temp = np.arange(2*viewOut_size+1, dtype=np.float32) + 1.0
    lorentz_row_temp = lorentz_row_temp - (viewOut_size+1)
    lorentz_col_temp = np.arange(2*viewOut_size+1, dtype=np.float32) + 1.0
    lorentz_col_temp = lorentz_col_temp - (viewOut_size+1)
    lorentz_col, lorentz_row = np.meshgrid(lorentz_col_temp, lorentz_row_temp)
#   Convert the frequency unit into lp/mm
    lorentz_col = lorentz_col/((cfg.scanner.detectorColSize/cfg.physics.FlatPanel_OSfactor)*len(lorentz_col_temp));
    lorentz_row = lorentz_row/((cfg.scanner.detectorRowSize/cfg.physics.FlatPanel_OSfactor)*len(lorentz_row_temp));
    H = cfg.scanner.FlatPanel_H
    lorentz_filt = 1 / (1 + H * (lorentz_col**2 + lorentz_row**2))
    
#   Perform convolution in frequency domain
    viewOut_freq = np.fft.fft2(viewOut, s=(2*viewOut_size+1, 2*viewOut_size+1))
    viewOut_filt_freq = np.fft.fftshift(viewOut_freq) * lorentz_filt
    viewOut_filt_spatial = abs(np.fft.ifft2(np.fft.ifftshift(viewOut_filt_freq)))
    viewOut_final = viewOut_filt_spatial[0:viewOut_size, 0:viewOut_size]
    
    viewOut_final_2 = viewOut_final.reshape(viewIn.shape)
    viewOut_final_2 = viewOut_final_2.astype('float32')
    
    return viewOut_final_2
