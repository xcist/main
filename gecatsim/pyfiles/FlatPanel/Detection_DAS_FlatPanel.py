# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import numpy as np
from scipy import signal

def Detection_DAS_FlatPanel(viewIn, viewId, cfg):
    '''
    - Performs image/frequency-domain down-sampling to model detector aliasing in the flat panel detector
    - Models the electronic gain conversion stage and the additive electronic noise

    Modified by: Junyuan Li, Department of Biomedical Engineering, Johns Hopkins University
    '''
    
    viewOut = viewIn*cfg.scanner.detectionGain
    
    # Additive electronic noise
    if cfg.physics.enableElectronicNoise:
        eNoise = np.float32(np.random.randn(viewIn.size)*cfg.sim.eNoise)
        viewOut += eNoise.reshape(viewIn.shape)
    
    # Model detector aliasing via image- or frequency-domain apparoach, given user's choice
    viewOut_size = int(cfg.scanner.detectorColCount*cfg.physics.FlatPanel_OSfactor)
    
    if cfg.physics.DetectorAliasMode == 0:
        # Modeling detector aliasing via image-domain approach (down-sampling process)
        viewOut = viewOut.reshape((viewOut_size, viewOut_size))
        SignalIntegral_mask = 1/(cfg.physics.FlatPanel_OSfactor**2) * np.ones((cfg.physics.FlatPanel_OSfactor, cfg.physics.FlatPanel_OSfactor))
        viewOut_filt = signal.convolve2d(viewOut, SignalIntegral_mask, mode='valid')
        # Extract subset of the presampling image
        viewOut_alias = viewOut_filt[0:viewOut_size:cfg.physics.FlatPanel_OSfactor, 0:viewOut_size:cfg.physics.FlatPanel_OSfactor]
    
    else:
        # Modeling detector aliasing via frequency-domain approach
        viewOut = viewOut.reshape((viewOut_size, viewOut_size))
        # Build the detector pixel MTF (2D sinc function) in frequency domain
        freq_row_temp1 = np.arange(-viewOut_size, viewOut_size+1.0, 1.0)
        freq_row_temp2 = freq_row_temp1.reshape((1, 2*viewOut_size+1))
        # Convert the frequency unit into cyc/mm
        freq_row = freq_row_temp2/((cfg.scanner.detectorColSize/cfg.physics.FlatPanel_OSfactor)*(2*viewOut_size+1))
        DetectorPixel_MTF_u = np.sinc(cfg.scanner.detectorColFillFraction*cfg.scanner.detectorColSize*freq_row)
        DetectorPixel_MTF_u_rep = np.tile(DetectorPixel_MTF_u, [2*viewOut_size+1,1])
        DetectorPixel_MTF_v_rep = np.transpose(DetectorPixel_MTF_u_rep)
        DetectorPixel_MTF = DetectorPixel_MTF_u_rep * DetectorPixel_MTF_v_rep
#        SignalIntegral_mask = 1/(cfg.physics.usfactor**2) * np.ones((cfg.physics.usfactor, cfg.physics.usfactor))
#        DetectorPixel_MTF = np.fft.fftshift(np.fft.fft2(SignalIntegral_mask, s=(2*viewOut_size+1, 2*viewOut_size+1)))
        
        viewOut_freq = np.fft.fft2(viewOut, s=(2*viewOut_size+1, 2*viewOut_size+1))
        viewOut_filt_freq = np.fft.fftshift(viewOut_freq) * DetectorPixel_MTF
        viewOut_filt_spatial = np.real(np.fft.ifft2(np.fft.ifftshift(viewOut_filt_freq)))
        viewOut_filt = viewOut_filt_spatial[1:(viewOut_size+1), 1:(viewOut_size+1)]
        
        # Downsampling in frequency domain
        viewOut_freq_noZeroPad = np.fft.fftshift(np.fft.fft2(viewOut_filt))
        freq_shift = np.arange(1, int(np.ceil(cfg.physics.FlatPanel_OSfactor)))
        # Superposition of shifted Fourier spectrum on x/y axis
        ScalingFactor = 1 / (cfg.physics.FlatPanel_OSfactor**2)
        viewOut_alias_freq = ScalingFactor * viewOut_freq_noZeroPad
    
        for Shift_axisXY in freq_shift:
            freq_shiftdist = int(Shift_axisXY * (1/cfg.physics.FlatPanel_OSfactor) * viewOut_size)
            spectrum_superimpose_1 = ScalingFactor * np.roll(viewOut_freq_noZeroPad, freq_shiftdist, axis=1)
            spectrum_superimpose_1[:, 0:freq_shiftdist] = 0.0
            spectrum_superimpose_2 = ScalingFactor * np.roll(viewOut_freq_noZeroPad, -freq_shiftdist, axis=1)
            spectrum_superimpose_2[:, -freq_shiftdist:] = 0.0
            spectrum_superimpose_3 = ScalingFactor * np.roll(viewOut_freq_noZeroPad, freq_shiftdist, axis=0)
            spectrum_superimpose_3[0:freq_shiftdist, :] = 0.0
            spectrum_superimpose_4 = ScalingFactor * np.roll(viewOut_freq_noZeroPad, -freq_shiftdist, axis=0)
            spectrum_superimpose_4[-freq_shiftdist:, :] = 0.0
            viewOut_alias_freq = viewOut_alias_freq + spectrum_superimpose_1 + spectrum_superimpose_2 + spectrum_superimpose_3 + spectrum_superimpose_4
            # Delete temporary variables
            del spectrum_superimpose_1, spectrum_superimpose_2, spectrum_superimpose_3, spectrum_superimpose_4, freq_shiftdist
            
        # Superposition of shifted Fourier spectrum in four quandrants
        for Shift_axis in freq_shift:
            freq_shiftdist = int(Shift_axis * (1/cfg.physics.FlatPanel_OSfactor) * viewOut_size)
            temp1 = ScalingFactor * np.roll(viewOut_freq_noZeroPad, freq_shiftdist, axis=1)
            spectrum_superimpose_1 = np.roll(temp1, freq_shiftdist, axis=0)
            spectrum_superimpose_1[0:freq_shiftdist, 0:freq_shiftdist] = 0.0
            temp2 = ScalingFactor * np.roll(viewOut_freq_noZeroPad, freq_shiftdist, axis=1)
            spectrum_superimpose_2 = np.roll(temp2, -freq_shiftdist, axis=0)
            spectrum_superimpose_2[-freq_shiftdist:, 0:freq_shiftdist] = 0.0
            temp3 = ScalingFactor * np.roll(viewOut_freq_noZeroPad, -freq_shiftdist, axis=1)
            spectrum_superimpose_3 = np.roll(temp3, freq_shiftdist, axis=0)
            spectrum_superimpose_3[0:freq_shiftdist, -freq_shiftdist:] = 0.0
            temp4 = ScalingFactor * np.roll(viewOut_freq_noZeroPad, -freq_shiftdist, axis=1)
            spectrum_superimpose_4 = np.roll(temp4, -freq_shiftdist, axis=0)
            spectrum_superimpose_4[-freq_shiftdist:, -freq_shiftdist:] = 0.0
            viewOut_alias_freq = viewOut_alias_freq + spectrum_superimpose_1 + spectrum_superimpose_2 + spectrum_superimpose_3 + spectrum_superimpose_4
            # Delete temporary variables
            del temp1, temp2, temp3, temp4, spectrum_superimpose_1, spectrum_superimpose_2, spectrum_superimpose_3, spectrum_superimpose_4, freq_shiftdist
            
        # Truncate at actual Nyquist frequency and perform inverse Fourier Transform
        crop_StartIdx = round(viewOut_size/2) - round(cfg.scanner.detectorColCount/2)
        crop_EndIdx = crop_StartIdx + cfg.scanner.detectorColCount
        viewOut_alias_crop = viewOut_alias_freq[crop_StartIdx:crop_EndIdx, crop_StartIdx:crop_EndIdx]
        viewOut_alias = np.real(np.fft.ifft2(np.fft.ifftshift(viewOut_alias_crop)))

    # Data type conversion
    viewOut_alias = viewOut_alias.astype('float32')
    
    return viewOut_alias.ravel()
