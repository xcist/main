# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

import numpy as np
import numpy.matlib
import copy
import ctypes
import os
import matplotlib.pyplot as plt
import gecatsim as xc
from gecatsim.pyfiles.CommonTools import *
from gecatsim.dose.pyfiles.calcDetectorFlux import calcDetectorFlux
from gecatsim.dose.pyfiles.img2vol import img2vol
from gecatsim.dose.pyfiles.C_DD3Dose import  C_DD3Dose
from gecatsim.dose.pyfiles.DoseConv import DoseConv

def load_C_lib():

    recon_lib = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib")

    # load C/C++ lib
    ll = ctypes.cdll.LoadLibrary
    if os.name == "nt":
        lib_file = "Dose_Recon_Library_Windows64.dll"
    else:
        lib_file = "Dose_Recon_Library_Linux64.so"
    clib = ll(os.path.join(recon_lib, lib_file))
    
    return clib

def load_catsim_lib():
    lib_path = my_path.paths["lib"]
    my_path.add_dir_to_path(lib_path)

    # load C/C++ lib
    ll = ctypes.cdll.LoadLibrary
    if os.name == "nt":
        libFile = "libcatsim64.dll"
    else:
        libFile = "libcatsim.so"
    clib = ll(os.path.join(lib_path, libFile))
    
    return clib


def catdoserecon(configfilename = None,cfg = None,adjust = None,preadjust = None,silent = None): 
    ## check arguments and load catsim lib, note this is not dose lib
    doselib = load_catsim_lib()
    
    ## view angles
    dviewangle = 2 * np.pi / cfg.protocol.viewsPerRotation
    viewangles = (np.arange(cfg.protocol.startViewId, cfg.protocol.startViewId+cfg.protocol.viewCount)) * dviewangle * cfg.protocol.rotationDirection
    
    ## detector geometry
    cfg = feval(cfg.scanner.detectorCallback,cfg)
    #cfg0 = copy.deepcopy(cfg)
    ## detector coordinates
    nrdetcols = cfg.scanner.detectorColCount
    nrdetrows = cfg.scanner.detectorRowCount
    #total_n_cells = nrdetcols*nrdetrows
    det_xyz = np.zeros((3, cfg.det.totalNumCells))
    #n_modules = int(cfg.scanner.detectorRowCount/cfg.scanner.detectorRowsPerMod)
    # TODO: what is this? do we need to do this
    #modtypes = [0]*n_modules # suppose we only have one type
    #n_cells = [cfg.scanner.detectorColsPerMod*cfg.scanner.detectorRowsPerMod] # cells per module
    for m in np.arange(cfg.det.nMod):
        type_ind = cfg.det.modTypes[m]
        n_cells = cfg.det.nCells #[type_ind] in xcist we only have one type
        cell_ind = np.arange(cfg.det.startIndices[m], cfg.det.startIndices[m]+n_cells)
        det_xyz[:,cell_ind] = np.matlib.repmat(cfg.det.modCoords[m],n_cells,1).T +\
                np.array([cfg.det.uvecs[m],cfg.det.vvecs[m]]).T @ cfg.det.cellCoords[np.arange(n_cells)].T
    det_xyz = np.reshape(det_xyz, (3, nrdetrows,nrdetcols), order='F')
    det_xyz = np.single(det_xyz)
    det_x = np.squeeze(det_xyz[0,int(np.ceil(nrdetrows / 2)),:])
    det_y = np.squeeze(det_xyz[1,int(np.ceil(nrdetrows / 2)),:])
    det_z = np.squeeze(det_xyz[2,:,int(np.ceil(nrdetcols / 2))])
    
    ## X-ray source position
    cfg = feval(cfg.scanner.focalspotCallback,cfg)
    # source coordinate
    src_xyz = np.single(np.mean(cfg.src.samples,0))
    
    ## ray angles
    cfg = feval(cfg.physics.rayAngleCallback, cfg)
    
    ## spectrum
    cfg.sim.subViewCount = 1
    cfg.sim.isOffsetScan = 0
    cfg = feval(cfg.protocol.spectrumCallback, cfg)
    #cfg.Evec = cfg.spec.Evec
    
    ## X-ray flux for all cells and views
    sino = np.single(np.zeros((cfg.physics.energyCount, cfg.det.totalNumCells, cfg.protocol.viewCount)))
    ViewTime = cfg.protocol.rotationTime/cfg.protocol.viewsPerRotation
    SubviewTime = ViewTime/cfg.physics.viewSampleCount
    # TODO: what is start time
    cfg.time = (cfg.protocol.startViewId - 0.5) * ViewTime + 0.5 * SubviewTime
    #cfg.time = cfg.start_time + (cfg.protocol.startViewId - 1.5) * ViewTime + 0.5 * SubviewTime
    
    # TODO: what is this for 
    cfg.dose.doselib = load_C_lib()
    ViewIndex = 0
    for ViewNumber in np.arange(cfg.protocol.startViewId,cfg.protocol.stopViewId+1):
        # Cat_Common_System_Setup
        detFlux = calcDetectorFlux(ViewIndex,ViewNumber,cfg)
        sino[:,:,ViewIndex] = detFlux.T
        ViewIndex = ViewIndex + 1
    
    # add the flux in the inactive area
    sino = sino * cfg.scanner.detectorColSize * cfg.scanner.detectorRowSize/cfg.det.activeArea
    ## other parameters
    nrviews = cfg.protocol.viewCount
    xoffset = cfg.dose.xoffset
    yoffset = cfg.dose.yoffset
    zoffset = cfg.dose.zoffset
    
    zshifts = np.arange(nrviews) * cfg.protocol.tableSpeed * cfg.protocol.rotationTime / cfg.protocol.viewsPerRotation
    
    n_voxel = cfg.dose.nVoxel
    nrcols = n_voxel
    nrrows = n_voxel
    nrplanes = cfg.recon.sliceCount
    
    voxel_xy_size = cfg.recon.fov/n_voxel
    voxel_z_size = cfg.recon.sliceThickness
    
    ## read recon image and convert to Mu volume
    #   vol:      Mu, cm^-1
    #   dens_vol: density, g/cm^3
    #   mass_vol: mass, g/voxel
    #   [nrcols, nrrows, nrplanes]
    vol, dens_vol, mass_vol = img2vol(cfg)
    # To get the Mu of each energy bin, we need its     material composition.
    # Currently consider only two kinds of material, water and bone, and simply
    # differentiate them with a threshold of Mu. -- Mingye
    if not hasattr(cfg, 'dose.doserecon_phantom_material') :
        ind_bone = vol > cfg.dose.waterThreshold
        ind_water = np.logical_and(vol<cfg.dose.waterThreshold, vol > cfg.dose.airThreshold)
        mask_highZ = np.single(ind_bone)
    
    # PROJECTION AND CONVOLUTION
    dosevol = np.single(np.zeros((nrcols,nrrows,nrplanes)))
    mydc = myDC()
    if len(cfg.spec.Evec.shape) < 1: cfg.spec.Evec = np.expand_dims(cfg.spec.Evec,0)
    for EnergyIndex in np.arange(cfg.physics.energyCount):
        print('Process Energy Bin # %d/%d\r' % (EnergyIndex+1, cfg.physics.energyCount))
        ee = [cfg.spec.Evec[EnergyIndex]]
        mu_water = xc.GetMu('water',ee)
        this_dosevol = np.single(np.zeros((nrcols,nrrows,nrplanes)))
        # flux at this energy
        this_sino = sino[EnergyIndex,:,:]
        this_sino = np.reshape(this_sino, [cfg.det.totalNumCells, cfg.protocol.viewCount], order="F")
        # Mu at this energy, Mu with Hui's BHC: water: 0.2 (63.95keV)
        this_vol = np.copy(vol)
        if not hasattr(cfg, 'dose.doserecon_phantom_material') :
            this_vol[ind_bone] = this_vol[ind_bone] * xc.GetMu('bone',ee) / cfg.dose.muBone
            this_vol[ind_water] = this_vol[ind_water] * xc.GetMu('water',ee) / cfg.dose.muWater
        else:
            this_vol[this_vol > 0.002] = this_vol((this_vol > 0.002)) * xc.GetMu(cfg.dose.doserecon_phantom_material,ee) / xc.GetMu(cfg.dose.doserecon_phantom_material,63.95)
        # dose tracing
        this_dosevol_tmp, this_sino_tmp = C_DD3Dose(cfg, src_xyz[0],src_xyz[1],src_xyz[2],\
                np.int32(nrdetcols),np.int32(nrdetrows),det_x,det_y,det_z,\
                np.single(xoffset),np.single(yoffset),np.single(zoffset),\
                np.single(viewangles),np.single(zshifts),\
                np.int32(nrviews),this_sino,np.int32(nrcols),np.int32(nrrows),np.int32(nrplanes),\
                this_vol,this_dosevol,np.single(voxel_xy_size),np.single(voxel_z_size))
        #         testfile=['test_',num2str(ee),'.mat'];
        #         load(testfile);
        #         save(testfile,'this_dosevol');
        # apply convolution
        # TODO should we also implement method 1
        if cfg.dose.doConvol == 1:
            #if hasattr(cfg.dose,'method') and cfg.method == 1:
            #    scipy.io.loadmat('dosereconkernel_mthd1.mat')

            #    x = dosereconkernel[ee]
            #    this_dosevol = doseconvol(cfg,this_dosevol,this_vol,mu_water,voxel_xy_size,x)
            #else:
            # DoseRecon 2.0
            this_dosevol_tmp = DoseConv(cfg, mydc, ee, this_dosevol_tmp, this_vol, voxel_xy_size,mask_highZ)
        # convert to deposited energy, in keV
        this_dosevol_tmp = this_dosevol_tmp * ee
        # convert to Dose, in mGy (Gy=J/kg)
        if hasattr(cfg.dose,'convertTomGy') and cfg.dose.convertTomGy == 1:
            this_dosevol_tmp = convert_to_mGy(this_dosevol_tmp,mass_vol)
        dosevol = dosevol + this_dosevol_tmp
        #breakpoint()
    
    ## write dosevol to raw file
    WriteDoseFiles(cfg, '', dosevol, 1)
    
    return dosevol
    
def WriteDoseFiles(cfg = None,ViewNumberString = None,DoseVolume = None,VerboseLocal = None): 
    if not hasattr(cfg.dose, 'doseFileName') :
        cfg.dose.doseFileName = cfg.results_basename
    
    if hasattr(cfg.dose,'convertTomGy') and cfg.dose.convertTomGy == 1:
        ext = '.mGy'
    else:
        ext = '.keV'
    
    dosefilename = cfg.dose.doseFileName+ViewNumberString+ext
    #outdose = np.copy(DoseVolume, ctype
    doseout = np.single(np.copy(np.transpose(DoseVolume, (2,1,0)), order='C'))
    xc.rawwrite(dosefilename, doseout)
    #Headline(np.array([sprintf('Output written to %s.',dosefilename)]),VerboseLocal)
    print('Output written to %s.'%dosefilename)

## convert to Dose, in mGy (Gy=J/kg)
# abs_vol is in keV
# mass_vol is in g/cm^3
def convert_to_mGy(abs_vol = None,mass_vol = None): 
    dose = abs_vol / mass_vol * (1000.0 * 1.6e-19 * 1000.0 * 1000.0)
    mask = np.logical_or(np.isnan(dose), np.isinf(dose))
    dose[mask] = 0
    return dose

# for persisten variables in 
class myDC():
    def __init__(self):
        self.data_Compt_dE_frac = None 
        self.data_E_scatter = None 
        self.dosereconkernel = None
        self.ee_CdE = None
        self.ee_mE = None
