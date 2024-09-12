# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.randpf import randpf
from gecatsim.pyfiles.CommonTools import *
from scipy.io import loadmat
import copy

def Detection_PC(cfg, viewId, subViewId):

    Evec = cfg.sim.Evec
    
    if viewId == cfg.sim.startViewId and subViewId == 0:
        ### detection efficiency
        # detector prefilter
        cfg.sim.Wvec = feval(cfg.physics.prefilterCallback, cfg)
        
        # detector absorption
        detectorMu = GetMu(cfg.scanner.detectorMaterial, Evec)
        detEff = 1-np.exp(-0.1*cfg.scanner.detectorDepth/np.cos(cfg.det.betas)*detectorMu)
        
        np.multiply(cfg.sim.Wvec, detEff, out=cfg.sim.Wvec)
        
        ### spectral response - primary
        cfg = get_spectral_response(cfg)
        
        ### crosstalk kernel
        cfg = get_xtalk_kernel(cfg)
        
    ### accumulate subviews (flux)
    if subViewId == 0:
        cfg.thisView = copy.copy(cfg.thisSubView)
    else:
        cfg.thisView += cfg.thisSubView
    
    ### for final subview
    if subViewId == cfg.sim.subViewCount-1:
        ### apply energy-dependent detection efficiency
        thisView = cfg.thisView*cfg.sim.Wvec
    
        ### apply spectral response to primary
        thisView = thisView @ cfg.sim.res_mat
    
        # scatter cross-talk
        # if cfg.physics.crosstalkCallback:
        #     xtalkSubView = feval(cfg.physics.crosstalkCallback, cfg.thisView, cfg)
        
        # detector bins
        thisView = rebin_energy(cfg, thisView)

        # quantum noise
        if cfg.sim.enableQuantumNoise:
            thisView = randpf(thisView)

        cfg.thisView = thisView
        
    return cfg

def get_spectral_response(cfg):
    Evec = cfg.sim.Evec
    kV = (Evec[1]-Evec[0])*len(Evec)
    Dvec = np.arange(0.5,kV)

    resp_fname = my_path.find("response_matrix", cfg.scanner.detectionResponseFilename, "")
    data = loadmat(resp_fname)
    Evec0 = data['Evec0'].ravel()
    Dvec0 = data['Dvec0'].ravel()
    res_mat0 = data['res_mat'][:,0:int(kV)]
    
    res_mat = np.zeros([len(Evec),len(Dvec)],dtype=np.float32)
    for indD in range(len(Dvec)):
        res_mat[:,indD] = np.interp(Evec,Evec0,res_mat0[:,indD])
    
    cfg.sim.Dvec = Dvec
    cfg.sim.res_mat = res_mat
    
    return cfg

def get_xtalk_kernel(cfg):

    return cfg

def rebin_energy(cfg, thisView):
    if hasattr(cfg.scanner,'detectorSumBins') and cfg.scanner.detectorSumBins:
        rebinView = thisView.sum(axis=1)
        return rebinView

    Dvec = cfg.sim.Dvec
    binThreshold = np.array(cfg.scanner.detectorBinThreshold)

    nPix = thisView.shape[0]
    nBin = len(binThreshold)-1

    # use overlap to calc the bin window weightings
    binWindow = np.zeros((len(Dvec),nBin),dtype=thisView.dtype)
    for i in range(nBin):
        tmpBin = np.zeros(nBin)
        tmpBin[i] = 1
        tmpWindow = overlap(None, tmpBin, Dvec, binThreshold)
        binWindow[:,i] = tmpWindow

    rebinView = thisView @ binWindow

    
    # rebinView1 = np.zeros((nPix,nBin),dtype=thisView.dtype)
    # for i in range(nPix):
    #     tmpSpec = thisView[i,:]
    #     rebinSpec = overlap(Dvec, tmpSpec, None, None, binThreshold)
    #     binWidth = np.diff(binThreshold)
    #     rebinSpec *= binWidth
    #     rebinView1[i,:] = rebinSpec

    # d = rebinView[2500,:]/rebinView1[2500,:]-1

    return rebinView
    
# if __name__ == "__main__":
#     pass
    
    
