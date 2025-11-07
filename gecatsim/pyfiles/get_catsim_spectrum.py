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
    spec = feval(cfg.protocol.spectrumCallback, cfg)
    if sim_step == 0:
        ss = spec.Ivec
        return ss, spec, cfg

    # Spectrum after bowtie & flat filter, before detector prefilter
    det0 = feval(cfg.scanner.detectorCallback, cfg)
    src0 = feval(cfg.scanner.focalspotCallback, cfg)
    cfg = Detector_RayAngles_2D(cfg)
    cfg.Evec = spec.Evec
    bowtie = feval(cfg.scanner.bowtieCallback, cfg)
    FiltrationTransVec = feval(cfg.scanner.filtrationCallback, cfg)
    spec = Combine_Spectrum_Bowtie_FlatFilter(cfg, bowtie, spec, FiltrationTransVec)
    if sim_step == 1:
        ss = spec.netIvec
        return ss, spec, cfg

    # Spectrum into detector active region (after detector prefilter)
    if hasattr(cfg.scanner, 'detectorPrefilterCallback'):
        Eff_prefilter = feval(cfg.scanner.detectorPrefilterCallback, cfg)
        spec.netIvec = spec.netIvec * Eff_prefilter
    if sim_step == 2:
        ss = spec.netIvec
        return ss, spec, cfg

    return None, spec, cfg