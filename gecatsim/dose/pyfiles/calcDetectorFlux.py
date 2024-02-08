# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

import numpy as np
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Gantry_Helical import Gantry_Helical
from gecatsim.dose.pyfiles.Combine_Spectrum_Bowtie_FlatFilter import Combine_Spectrum_Bowtie_FlatFilter

def calcDetectorFlux(ViewIndex, ViewNumber, cfg):
    # source
    if ViewIndex == 0 or cfg.physics.recalcSrc:
        cfg = feval(cfg.scanner.focalspotCallback,cfg)
        # Note the initial source position.
        #cfg.src0 = src0
        # TODO it seems that we do not have weights in setsource
        if cfg.phantom.callback == "Phantom_Analytic":
            cfg.callback_setSource = 'C_Source_Analytic_Set';
        elif cfg.phantom.callback == "Phantom_NCAT":
            cfg.callback_setSource = 'C_Source_NCAT_Set';
        elif cfg.phantom.callback == "Phantom_Voxelized":
            cfg.callback_setSource = 'C_Source_Voxelized_Set';
        #feval(cfg.callback_setSource,cfg, cfg.src.weights,cfg.src.nSamples)
        cfg.src.weights = feval(cfg.callback_setSource,cfg, cfg.src.weights,cfg.src.nSamples)

    # spectrum
    # Setup spectra for Variable kVp & mA simulation, Mingye
    if ViewIndex == 0 and hasattr(cfg,'variable_kV_mA_sim') and cfg.variable_kV_mA_sim == 1:
        # read spectra for all kVps
        kVs = np.arange(np.min(cfg.kV_profile),np.max(cfg.kV_profile)+1)
        spec_all_kV = cell(np.asarray(kVs).size,1)
        for kV_id in np.arange(np.asarray(kVs).size):
            cfg.spectrum_filename = sprintf_filepath(cfg.spectrum_base_filename,kVs(kV_id))
            cfg.number_Ebins = kVs(kV_id) / cfg.Ebin_sampling_factor
            spec_all_kV[kV_id],cfg_tmp = feval(cfg.callback_spectrum,cfg)
        # set E-bins and Evec to max kV for all views
        cfg = cfg_tmp
        spec = spec_all_kV[end()]
        # fullfil the high empty Ebins with zeros
        Ivec_all_kV = cell(np.asarray(kVs).size,1)
        for kV_id in np.arange(1,np.asarray(kVs).size+1).reshape(-1):
            dim = spec_all_kV[kV_id].Ivec.shape
            dim[1] = spec.number_Ebins - dim(1)
            Ivec_all_kV[kV_id] = np.array([[spec_all_kV[kV_id].Ivec],[np.zeros((dim,dim))]])
    
    if ViewIndex == 0 or cfg.physics.recalcSpec:
        # Variable kVp & mA simulation, Mingye
        if hasattr(cfg,'variable_kV_mA_sim') and cfg.variable_kV_mA_sim == 1:
            kV_id = find(kVs == cfg.kV_profile(ViewNumber))
            #fprintf('kV_id #d ',kV_id);
            mA_factor = 1
            if hasattr(cfg,'mA_profile') and np.asarray(cfg.mA_profile).size == cfg.total_n_views:
                mA_factor = cfg.mA_profile(ViewNumber) / cfg.mA
            spec.Ivec = Ivec_all_kV[kV_id] * mA_factor
            cfg.Ivec = spec.Ivec
        else:
            # constant kVp simulation
            cfg = feval(cfg.protocol.spectrumCallback, cfg)
    
    # bowtie, no such thing
    #if ViewIndex == 0:# or cfg.recompute_bowtie:
    #    bowtie = feval(cfg.callback_bowtie,cfg)
    
    # flat filters
    if ViewIndex == 0 or cfg.physics.recalcFilt:
        if cfg.protocol.filterCallback =='Filter_Flat_Dynamic':
            cfg = feval(cfg.protocol.filterCallback,cfg,ViewIndex)
        else:
            cfg = feval(cfg.protocol.filterCallback,cfg)
    
    # Combine SPECTRUM, BOWTIE, and FLAT FILTER
    if ViewIndex == 0 or cfg.physics.recalcSpec or cfg.physics.recalcFilt:
        cfg = Combine_Spectrum_Bowtie_FlatFilter(cfg)
    
    # GANTRY transforms
    # Compute & apply gantry transforms (relative to the initial source and detector positions).
    cfg = Gantry_Helical(cfg,ViewIndex)
    # DETECTOR FLUX
    if (ViewIndex == 0) or cfg.physics.recalcFlux:
        if len(cfg.physics.fluxCallback)!=0:
            #DetectorFlux = feval(cfg.callback_flux,spec,src,det,cfg)
            cfg = feval(cfg.physics.fluxCallback, cfg)
        else:
            # default detector signal if correct noise model is not required
            cfg.detFlux = 100000.0 * np.ones((spec.number_Ebins,det.total_n_cells))
    
    # mA MODulation
    if hasattr(cfg, 'callback_mA_modulation') and len(cfg.callback_mA_modulation)!=0:
        mA_modulation = feval(cfg.callback_mA_modulation,cfg, ViewNumber)
    else:
        mA_modulation = 1

    return  cfg.detFlux
