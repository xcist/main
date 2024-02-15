import os
import copy
import numpy as np
from gecatsim.pyfiles.CommonTools import rawread

def Prep_BHC_Accurate(cfg, prep):
    print("Applying Beam Hardening Correction (ACCURATE BHC)...", end='')
    
    if hasattr(cfg.physics, "BHC_vec_fname"):
        if os.path.isfile(cfg.physics.BHC_vec_fname):
            poly_coef = np.load(cfg.physics.BHC_vec_fname)
        else:
            poly_coef = gen_BHC_vec(cfg)
    else:
        poly_coef = gen_BHC_vec(cfg)

    num_poly_coefs = cfg.physics.BHC_poly_order+1
    for viewId in range(cfg.protocol.viewCount):
        view_out = poly_coef[:, 0]
        for i in range(1, num_poly_coefs):
            view_out = view_out*prep[viewId] + poly_coef[:,i]

        prep[viewId] = view_out

    print("done.\n")
    return prep

def gen_BHC_vec(cfg):

    from gecatsim.pyfiles.CatSim import CatSim
    ct = CatSim()
    ct.cfg_to_self(cfg)

    poly_order = ct.physics.BHC_poly_order
    num_poly_coefs = poly_order+1

    if not hasattr(ct.physics, "BHC_material") or ct.physics.BHC_material == "":
        mt = "water"
    else:
        mt = ct.physics.BHC_material

    max_length_mm = ct.physics.BHC_max_length_mm
    length_step_mm = ct.physics.BHC_length_step_mm
    num_length = int(max_length_mm/length_step_mm)

    ndet = ct.scanner.detectorColCount*ct.scanner.detectorRowCount
    sig = np.zeros((num_length, ndet))

    mtViews = rawread(ct.resultsName+".air", [ndet], 'float')
    I0 = mtViews

    orginal_prefilter = copy.copy(ct.scanner.detectorPrefilter)
    for j in range(num_length):
        thick = (1+j)*length_step_mm

        ct.scanner.detectorPrefilter = orginal_prefilter + [mt, thick]
        ct.air_scan(doPrint=False)
        rawViews = rawread(ct.resultsName+".air", [ndet], 'float')
        sig[j] = -np.log(rawViews/I0)

    ct.scanner.detectorPrefilter = orginal_prefilter
    ct.air_scan()

    dessig = ct.physics.EffectiveMu*length_step_mm/10*(1+np.arange(num_length))
    poly_coef = np.zeros((ndet, num_poly_coefs))
    for n in range(ndet):
        poly_coef[n] = np.polyfit(sig[:,n], dessig, poly_order)

    return poly_coef
