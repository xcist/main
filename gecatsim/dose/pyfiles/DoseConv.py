# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

import numpy as np
from scipy import interpolate,io
import gecatsim as xc
from gecatsim.pyfiles.CommonTools import *
from gecatsim.dose.pyfiles.GetMuByProcess import GetMuByProcess
from gecatsim.dose.pyfiles.get_voxel_stopping_power import get_voxel_stopping_power
from gecatsim.dose.pyfiles.doseconvol import doseconvol

def DoseConv(cfg = None, mydc=None, ee = None,dosevol_int = None,mu_vol = None,voxel_xy_size = None,mask_highZ = None): 
    mask_lowZ = 1 - mask_highZ
    
    if mydc.data_Compt_dE_frac is None or mydc.data_E_scatter is None or mydc.dosereconkernel is None:
        # Geant4 results: 'Compton_ELoss_factor',ee=10:10:160
        scat_dir = my_path.find("dose_data", 'scatter_E_loss_fraction.mat', '')
        scat_data = io.loadmat(scat_dir)
        mydc.data_Compt_dE_frac = np.squeeze(np.vstack(([[0]],scat_data['Compton_ELoss_factor'])))
        mydc.ee_CdE = np.array([0]+list(np.arange(10,160+10,10)))
        # mean scatter energy of monos, weighted by scatter spectrum, of multi-material
        # mean_E, mspr, Geant4 results,ee=10:10:160
        # load weighted_scatter_energy_by_spec.mat;
        tmp_dir2 = my_path.find("dose_data", 'weighted_scatter_energy_by_spec.mat', '')
        tmpdata2 = io.loadmat(tmp_dir2)
        mydc.data_E_scatter = np.vstack(([0], np.mean(tmpdata2['mean_E'],1,keepdims=True)))
        mydc.ee_mE = np.array([0]+list(np.arange(10,160+10,10)))
        # convolution kernel
        #load dosereconkernel_mthd4.mat; # dosereconkernel

        tmp_dir3 = my_path.find("dose_data", 'dosereconkernel_mthd4.mat', '')
        tmpdata3 = io.loadmat(tmp_dir3)
        mydc.dosereconkernel = tmpdata3['dosereconkernel']
    
    ## calc scattered energy
    # fraction of scattered energy to interaction energy
    f_interp = interpolate.interp1d(mydc.ee_CdE, mydc.data_Compt_dE_frac)
    Compt_dE_frac = f_interp(ee)
    frac_water = (xc.GetMu('water',ee) - GetMuByProcess(cfg,'water',ee,'PE') - GetMuByProcess(cfg,'water',ee,'C') * Compt_dE_frac) / xc.GetMu('water',ee)
    frac_bone = (xc.GetMu('bone',ee) - GetMuByProcess(cfg,'bone',ee,'PE') - GetMuByProcess(cfg,'bone',ee,'C') *Compt_dE_frac) / xc.GetMu('bone',ee)
    frac_vol = frac_water * mask_lowZ + frac_bone * mask_highZ
    # apply mu_frac to dosevol (interaction energy volume)
    dosevol_scatter = np.multiply(dosevol_int,frac_vol)
    dosevol_1st_abs = dosevol_int - dosevol_scatter
    ## convolution on dosevol_scatter
    # stopping power to scatters
    f_interp = interpolate.interp1d(mydc.ee_mE, np.squeeze(mydc.data_E_scatter))
    E_scatter = f_interp(ee)
    stp_water = get_voxel_stopping_power(cfg,'water',E_scatter,voxel_xy_size,1)
    stp_bone = get_voxel_stopping_power(cfg, 'bone',E_scatter,voxel_xy_size,1)
    # absorbing vol is determined by material stopping power and density
    dens_vol = mu_vol / (mask_lowZ * xc.GetMu('water',ee) + mask_highZ * xc.GetMu('bone',ee))
    abs_vol = np.multiply(dens_vol,(mask_lowZ * stp_water + mask_highZ * stp_bone))
    # convolution
    x = mydc.dosereconkernel[np.int32(np.ceil(ee))-1]
    x = x[0,0]
    #breakpoint()
    dosevol_2nd_abs = doseconvol(cfg,dosevol_scatter,abs_vol,stp_water,voxel_xy_size,x)
    ## final absorption
    dosevol_final = dosevol_1st_abs + dosevol_2nd_abs
    return dosevol_final

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
    import matplotlib.pyplot as plt
    from gecatsim.pyfiles.CommonTools import *
    
    matin = sio.loadmat("unittest/doseconpar_in.mat")
    cfg = CFG()
    mydc = myDC()
    
    pyout = DoseConv(cfg, mydc, matin['ee'][0], matin['dosevol_int'], matin['mu_vol'], matin['voxel_xy_size'][0,0], matin['mask_highZ'])

    matout = sio.loadmat("unittest/doseconpar_out.mat")
    matout = matout['dosevol_final']

    try:
        assert np.allclose(pyout, matout, atol=1.E-8), 'dosevol'
    except AssertionError as err:
        print(err)
