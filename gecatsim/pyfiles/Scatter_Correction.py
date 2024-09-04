# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy.matlib as nm
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Scatter_ConvolutionModel import get_scatter_kernel
from tqdm import tqdm

# This is a simplified kernel based scatter correction algorithm.
def Scatter_Correction(cfg, airscan, offsetScan, phantomScan):
    print("Applying Scatter Correction...", end='')

    ###--------- Get scatter kernel
    if cfg.physics.scatterKernelCallback:
        cfg.scatter_kernel = feval(cfg.physics.scatterKernelCallback, cfg)
    else:
        cfg.scatter_kernel = get_scatter_kernel()
            
    ###--------- -log
    #if cfg.protocol.airViewCount==1:
    #    airscan = nm.repmat(airscan, cfg.protocol.viewCount, 1)
    #if cfg.protocol.offsetViewCount==1:
    #    offsetScan = nm.repmat(offsetScan, cfg.protocol.viewCount, 1)
    prep = (phantomScan-offsetScan)/(airscan-offsetScan)
    smallValue = 1.E-10
    prep[prep<smallValue] = smallValue
    prep = -np.log(prep)
    prep[prep<smallValue] = smallValue

    for viewId in range(cfg.protocol.viewCount):
        if not hasattr(cfg.physics, "scatterCorrectionScaleFactor"):
            cfg.physics.scatterCorrectionScaleFactor = 1
        sc_preConv = phantomScan[viewId,:]*np.power(prep[viewId,:],0.9)*0.0268*cfg.physics.scatterCorrectionScaleFactor
        sc_preConv = sc_preConv.reshape(cfg.scanner.detectorRowCount, cfg.scanner.detectorColCount)
        sc_conv = conv2(sc_preConv, cfg.scatter_kernel, 'same')
        sc_conv = sc_conv.ravel()
        phantomScan[viewId,:] -= sc_conv
        
        # import matplotlib.pyplot as plt
        # sc_conv = sc_conv.reshape(cfg.scanner.detectorRowCount, cfg.scanner.detectorColCount)
        # plt.plot(sc_conv[32, :])
        # plt.show()
    
    print("done.\n")
    
    ###--------- save corrected scan
    if hasattr(cfg.physics,"scatterCorrectionSaveView") and cfg.physics.scatterCorrectionSaveView:
        rawwrite(cfg.resultsName+'_SC.scan', phantomScan)
    
    
    return airscan, offsetScan, phantomScan

    
