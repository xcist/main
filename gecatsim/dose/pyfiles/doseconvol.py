# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

import numpy as np
from gecatsim.dose.pyfiles.brconvol_matlab import brconvol_matlab

# input : dosevol, vol, mu_water, voxel_xy_size, x(kernel) 
def doseconvol(cfg = None,dosevol = None,vol = None,mu_water = None,voxel_xy_size = None,x = None): 
    f0 = dosevol.flatten()
    f1 = np.multiply(brconvol_matlab(dosevol,0.85 / voxel_xy_size,-1),(vol / mu_water))
    f1 = f1.flatten()
    f4 = np.multiply(brconvol_matlab(dosevol,3.5 / voxel_xy_size,-1),(vol / mu_water))
    f4 = f4.flatten()
    f8 = np.multiply(brconvol_matlab(dosevol,8.0 / voxel_xy_size,-1),(vol / mu_water))
    f8 = f8.flatten()
    f16 = np.multiply(brconvol_matlab(dosevol,16.0 / voxel_xy_size,-1),(vol / mu_water))
    f16 = f16.flatten()
    f32 = np.multiply(brconvol_matlab(dosevol,32.0 / voxel_xy_size,-1),(vol / mu_water))
    f32 = f32.flatten()
    f64 = np.multiply(brconvol_matlab(dosevol,64.0 / voxel_xy_size,-1),(vol / mu_water))
    f64 = f64.flatten()
    f128 = np.multiply(brconvol_matlab(dosevol,128.0 / voxel_xy_size,-1),(vol / mu_water))
    f128 = f128.flatten()
    #[x,resnorm]=lsqlin(double([f0,f1,f4,f8,f16,f32, f64,f128]),double(Edep(:)),[],[],[1,1,1,1,1,1,1,1],[1.0],[0,0,0,0,0,0,0,0]',[]);
    Efit = np.array([f0,f1,f4,f8,f16,f32,f64,f128]).T@x
    dosevol_rt = np.reshape(Efit, dosevol.shape)
    return dosevol_rt

if __name__ == '__main__':
    # for persisten variables in 
    class myDC():
        def __init__(self):
            self.data_Compt_dE_frac = None 
            self.data_E_scatter = None 
            self.dosereconkernel = None
            self.ee_CdE = None
            self.ee_mE = None
            
    from scipy import io as sio
    
    matin = sio.loadmat("unittest/doseconv_in.mat")
    
    mydc = myDC()
    pyout = doseconvol(mydc, matin['dosevol'], matin['vol'], matin['mu_water'], matin['voxel_xy_size'], matin['x'])

    matout = sio.loadmat("unittest/doseconv_out.mat")
    matout = matout['dosevol']

    try:
        assert np.allclose(pyout, matout, atol=1.E-8), 'dosevol'
    except AssertionError as err:
        print(err)

