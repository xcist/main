# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

from ctypes import *
from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import *

def Phantom_Voxelized(cfg):
    
    ###----------- phantom file (vp file is json format)
    cfg.phantom.filename = my_path.find('phantom', cfg.phantom.filename, '')
    filepath, vpFilename = os.path.split(cfg.phantom.filename)
    with open(cfg.phantom.filename) as fin:
        vp = json.load(fin)
        
    ###----------- pass material Mu to C
    materialList = vp['mat_name']
    numberOfMaterials = vp['n_materials']
    cfg.phantom.numberOfMaterials = numberOfMaterials
    set_material(cfg, materialList)
    
    ###----------- pass volume info to C
    # import matplotlib.pyplot as plt # Only needed for plots of phantom (below)

    for i in range(numberOfMaterials):
        dims = np.array([vp['cols'][i], vp['rows'][i], vp['slices'][i]], dtype=np.int32)
        offsets = np.array([vp['x_offset'][i], vp['y_offset'][i], vp['z_offset'][i]], dtype=np.single)
        voxelsize = np.array([vp['x_size'][i], vp['y_size'][i], vp['z_size'][i]], dtype=np.single)*cfg.phantom.scale
        voxelsizeXY = voxelsize[0]
        voxelsizeZ  = voxelsize[2]
        
        volumePathName = cfg.phantom.filename.replace(vpFilename, vp['volumefractionmap_filename'][i])
        volumeData = rawread(volumePathName, [], vp['volumefractionmap_datatype'][i])
        # if density_scale is set
        if 'density_scale' in vp:
            volumeData *= vp['density_scale'][i]
        
        '''
        Offsets in .vp file indicate the index coordinate of the scanner origin (isocenter) relative to the edge of the phantom.
        DD Proj expects an offset (in mm) between the phantom voxel grid center and the scanner isocenter
        Here we fix the offsets to be what DD Proj expects
        '''
        offsets = ((np.single(dims)+1)/2-offsets)*voxelsize + np.array(cfg.phantom.centerOffset, dtype=np.single)

        ###------- prepair data for setting volume in C
        volumeDims = np.append(dims, np.int32(numberOfMaterials))
        xyMask = np.ones((dims[1]*2, dims[0]), dtype=np.uint8)
        
        # The order of bin data of CatSim phantom file is x -> y -> z
        # when python reads it, the matrix dim becomes [z, y, x]
        # And the C func wants data order: z -> x -> y
        # Note: numpy array starts from the last dim, and when transposed, .ravel is needed, or the order will not change.
        
        volumeData = volumeData.reshape(dims[2], dims[1], dims[0])
        volumeData = np.transpose(volumeData, (1, 2, 0)).ravel()
        volumeData = volumeData.astype(np.float32)

        # # Plot slice 0 for each material
        # plt.figure(i)
        # plt.imshow(volumeData.reshape(dims[2], dims[1], dims[0])[0], cmap='gray')
        # plt.title(materialList[i])

        materialIndex = i+1        
        set_voxelized_volume(cfg, volumeData, volumeDims, offsets, voxelsize, xyMask, materialIndex, numberOfMaterials)
    
    # plt.pause(1)
    # print('*******************************************')
    # print('* Press Enter to close plots and continue *')
    # input('*******************************************')
    # plt.close('all')

    return cfg

def set_material(cfg, materialList):
    Evec = cfg.spec.Evec
    nMat = len(materialList)
    Mus = np.zeros([Evec.size, nMat], dtype=np.single)
    for i in range(nMat):
        Mus[:, i] = GetMu(materialList[i], Evec)/10 # cm^-1 --> mm^-1
        
    # the C func wants data order: materialIndex -> Ebin, so the dim is [Ebin, materialIndex]
    fun = cfg.clib.set_material_info_vox
    fun.argtypes = [c_int, c_int, ndpointer(c_float)]
    fun.restype = None
    fun(nMat, Evec.size, Mus)

def set_voxelized_volume(cfg, volumeData, volumeDims, offsets, voxelsize, xyMask, materialIndex, numberOfMaterials):
    if hasattr(cfg.phantom, "useUInt16") and cfg.phantom.useUInt16:
        fun = cfg.clib.set_phantom_info_vox_uint16
        fun.argtypes = [POINTER(c_int), ndpointer(c_ushort), ndpointer(c_int), \
            c_float, c_float, c_float, c_float, c_float, ndpointer(c_ubyte), c_int, c_int]

        # convert data to uint16
        volumeData = np.ushort(np.round(volumeData*10000))
    else:
        fun = cfg.clib.set_phantom_info_vox
        fun.argtypes = [POINTER(c_int), ndpointer(c_float), ndpointer(c_int), \
            c_float, c_float, c_float, c_float, c_float, ndpointer(c_ubyte), c_int, c_int]
    fun.restype = None
    
    Status = [0]
    Status = (c_int*1)(*Status)
    fun(Status, volumeData, volumeDims, \
        offsets[0], offsets[1], offsets[2], voxelsize[0], voxelsize[2], xyMask, materialIndex, numberOfMaterials)
    
    if Status[0] == -1:
        print("*** Error %d in set_phantom_info_vox: memory allocation error!" % Status[0])
    elif Status[0] == -2:
        print("*** Error %d in set_phantom_info_vox: Not have enough system memory for voxelized phantom!" % Status[0])


# if __name__ == "__main__":
#     cfg = source_cfg("./cfg/default.cfg")
#     cfg.phantom.filename = 'BIG.json'
#
#     cfg = feval(cfg.scanner.detectorCallback, cfg)
#     cfg = feval(cfg.scanner.focalspotCallback, cfg)
#     cfg = feval(cfg.protocol.spectrumCallback, cfg)
#
#     cfg = Phantom_Voxelized(cfg)
