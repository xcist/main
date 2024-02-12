# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

import numpy as np
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Gantry_Helical import Gantry_Helical

def calcDetectorFlux(viewId, ViewNumber, cfg):
    # source
    if viewId == 0 or cfg.physics.recalcSrc:
        cfg = feval(cfg.scanner.focalspotCallback, cfg)

    # spectrum
    # Setup spectra for Variable kVp & mA simulation, Mingye
    if viewId == 0 and hasattr(cfg,'variable_kV_mA_sim') and cfg.variable_kV_mA_sim == 1:
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
    
    if viewId == 0 or cfg.physics.recalcSpec:
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
    
    # filters (bowtie and flat)
    if viewId == 0 or cfg.physics.recalcFilt:
        cfg = feval(cfg.protocol.filterCallback, cfg)
    
    # flux
    if viewId == 0 or cfg.physics.recalcFlux:
        cfg = feval(cfg.physics.fluxCallback, cfg)

    return  cfg.detFlux
