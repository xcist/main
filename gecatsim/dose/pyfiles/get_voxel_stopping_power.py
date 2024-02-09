# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

import numpy as np
from scipy import interpolate,io
import gecatsim as xc
from gecatsim.pyfiles.CommonTools import *
from gecatsim.dose.pyfiles.GetMuByProcess import GetMuByProcess
    
def get_voxel_stopping_power(cfg, material = None,ee = None,voxel_size = None,dens = None): 
    # Geant4 results: 'Compton_ELoss_factor',ee=10:10:160
    #load scatter_E_loss_fraction.mat;
    scat_dir = my_path.find("dose_data", 'scatter_E_loss_fraction.mat', '')
    scat_data = io.loadmat(scat_dir)
    data_Compt_dE_frac = np.squeeze(scat_data['Compton_ELoss_factor'])
    # interpolation to Evec
    #breakpoint()
    #ee_vec = np.vstack(([0],np.arange(10,160+10,10),[200]))
    ee_vec = np.array([0]+list(np.arange(10,160+10,10))+[200])
    #cf_vec = np.vstyack(([[0],[data_Compt_dE_frac],[0]]))
    cf_vec = np.array([0]+list(data_Compt_dE_frac)+[0])
    
    cf_vec[-1] = (200 - 160) / 20 * (cf_vec[-2] - cf_vec[-4]) + cf_vec[-2]
    # shortcut
    d = voxel_size * dens
    # 1st interaction
    mu_pe = GetMuByProcess(cfg,material,ee,'PE')
    mu_c = GetMuByProcess(cfg,material,ee,'C')
    mu_r = GetMuByProcess(cfg,material,ee,'R')
    mu = np.array(xc.GetMu(material,ee))
    f_interp = interpolate.interp1d(ee_vec,cf_vec)
    cf = f_interp(ee)
   
    r1 = np.multiply((1 - np.exp(- mu * d)),(mu_pe + np.multiply(cf,mu_c))) / mu
    ratio = r1
    return ratio
