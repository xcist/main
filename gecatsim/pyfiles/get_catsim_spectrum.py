# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import os
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Detector_RayAngles_2D import Detector_RayAngles_2D
from gecatsim.dose.pyfiles.Combine_Spectrum_Bowtie_FlatFilter import Combine_Spectrum_Bowtie_FlatFilter

def get_catsim_spectrum(cfg, sim_step=0):
    """
    cfg: the cfg that gets into catsim
    sim_step = 0: tube spectrum
               1: spectrum after bowtie and flat filters
               2: spectrum after detector prefilters
    ss: the spectrum, dim: [Ebin cfg.total_n_cells]
        unit: photons/view/mm^2 per energy bin at 1m
    """

    # Spectrum before bowtie & flat filter
    spec = cfg.callback_spectrum(cfg)
    if sim_step == 0:
        ss = spec.Ivec
        return ss, spec, cfg

    # Spectrum after bowtie & flat filter, before det prefilter
    det0 = cfg.callback_detector(cfg)
    src0 = cfg.callback_source(cfg)
    cfg = Detector_RayAngles_2D(cfg, det0, src0)
    cfg.Evec = spec.Evec
    bowtie = cfg.callback_bowtie(cfg)
    FiltrationTransVec = cfg.callback_filtration(cfg)
    spec = Combine_Spectrum_Bowtie_FlatFilter(cfg, bowtie, spec, FiltrationTransVec)
    if sim_step == 1:
        ss = spec.netIvec
        return ss, spec, cfg

    # Spectrum into detector active region (after detector prefilter)
    if hasattr(cfg, 'callback_detector_prefilter'):
        Eff_prefilter = cfg.callback_detector_prefilter(cfg)
        spec.netIvec = spec.netIvec * Eff_prefilter
    if sim_step == 2:
        ss = spec.netIvec
        return ss, spec, cfg

    return None, spec, cfg
