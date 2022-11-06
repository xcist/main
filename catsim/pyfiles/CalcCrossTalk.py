# Copyright 2022, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#import sys
#import numpy.matlib as nm
#from scipy import interpolate, io
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
#from catsim.pyfiles.CommonTools import *

def getConvKernel(row_crosstalk, col_crosstalk):
    row_ker = np.array([row_crosstalk, 1.-2.*row_crosstalk, row_crosstalk])
    col_ker = np.array([col_crosstalk, 1.-2.*col_crosstalk, col_crosstalk])
    return col_ker[:,None]*row_ker

# for xray xtalk, it is 2d
def CalcCrossTalk(thisView, cfg):

    # thisView is of dim: (row*channel)[xenergyCount]
    conv_kernel = getConvKernel(cfg.physics.row_crosstalk, cfg.physics.col_crosstalk)
    if conv_kernel[1,1] == 1.: return thisView
    
    outView = np.reshape(np.copy(thisView),
                [cfg.scanner.detectorColCount, cfg.scanner.detectorRowCount, cfg.physics.energyCount])
                #[cfg.scanner.detectorRowCount, cfg.scanner.detectorColCount, cfg.physics.energyCount])
    for i in range(outView.shape[2]):
        outView[:,:,i] = signal.convolve2d(outView[:,:,i], conv_kernel, mode='same', boundary='fill', fillvalue=0.)
    outView = outView.astype(np.single)

    #breakpoint()
    return np.reshape(outView, thisView.shape)

if __name__ == "__main__":
    import scipy.io
    import catsim as xc
    mat_in = scipy.io.loadmat('view_in_reshape.mat')['view_out']
    mat_out = scipy.io.loadmat('view_out_reshape.mat')['view_out']
    ct = xc.CatSim("all_in_one")
    cfg = ct.get_current_cfg()
    py_out = CalcCrossTalk(mat_in, cfg)
    diff = py_out-mat_out
    plt.imshow(diff[:,:32]);plt.colorbar();plt.show()
