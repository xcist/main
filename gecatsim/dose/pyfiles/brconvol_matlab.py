# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

import numpy as np
from scipy.signal import convolve, convolve2d
    
def brconvol_matlab(data = None,sigma = None,dim = None): 
    data = np.copy(data)
    if sigma != 0:
        # determine convolution directions
        if dim is None: dim = -1
        do_row = (dim == -1 or dim==255 or dim == 0)
        do_col = (dim == -1 or dim==255 or dim == 1)
        do_plane = (dim == -1 or dim==255 or dim == 2)
        kernelsize = np.max(data.shape) * 2 + 1
        kernel = np.exp(- (np.arange(kernelsize) - np.max(data.shape)) ** 2 / (2 * sigma ** 2))
        kernel = kernel.flatten()
        kernel /= np.sum(kernel)
        mid = np.max(data.shape)
        # convolve in row direction
        if (do_row):
            half = data.shape[1-1]
            k0 = kernel[np.arange(mid - half,mid + half+1)]
            #k0 = k0 / np.sum(k0)
            #for i in range(data.shape[2]):
            #    for j in range(data.shape[0]):
            #        data[j,:,i] = convolve(data[j,:,i], k0,'same')
            data = convolve(data, k0[None,:,None],'same')
        # convolve in col direction
        if (do_col):
            half = data.shape[2-1]
            k0 = kernel[np.arange(mid - half,mid + half+1)]
            #k0 = k0 / np.sum(k0)
            #data = np.convolve(data,np.transpose(k0),'same')
            #for i in range(data.shape[2]):
            #    for j in range(data.shape[1]):
                    #TODO: check this is the col ops
            #        data[:,j,i] = convolve(data[:,j,i], k0,'same')
            data = convolve(data, k0[:,None,None],'same')
        # convolve in plane direction
        if (do_plane):
            half = data.shape[3-1]
            k0 = kernel[np.arange(mid - half,mid + half+1)]
            #k0 = k0 / np.sum(k0)
            #data = np.convolve(data,permute(k0,np.array([1,2,0])),'same')
            if data.shape[2]==1:
                # k0 *= 5  # dose.method==1
                k0 *= 2  # DoseRecon 2.0
            #for i in range(data.shape[0]):
            #    for j in range(data.shape[1]):
            #        data[i,j,:] = convolve(data[i,j,:], k0, 'same')
            data = convolve(data, k0[None, None, :], 'same')
    
    return data

if __name__ == '__main__':
    from scipy import io as sio
    mat_in = sio.loadmat("unittest/brconv_in.mat")

    py_out = brconvol_matlab(mat_in['data'], mat_in['sigma'], mat_in['dim']-1)

    mat_out = sio.loadmat("unittest/brconv_out.mat")
    mat_out = mat_out['data']

    try:
        assert np.allclose(py_out, mat_out, atol=1.E-8), 'dosevol'
    except AssertionError as err:
        print(err)
