# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

from ctypes import *
from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import *


def Phantom_NCAT(cfg):
    ###----------- pass material Mu and volume to C
    nMat = set_material(cfg)
    set_volume(cfg, nMat)
    
    ###----------- pass detector and source info to C
    set_detector(cfg)
    set_source(cfg)

    return cfg


def set_material(cfg):
    materialList = [
        'ncat_water',
        'ncat_muscle',
        'ncat_lung',
        'ncat_dry_spine',
        'ncat_dry_rib',
        'ncat_adipose',
        'ncat_blood',
        'ncat_heart',
        'ncat_kidney',
        'ncat_liver',
        'ncat_lymph',
        'ncat_pancreas',
        'ncat_intestine',
        'ncat_skull',
        'ncat_cartilage',
        'ncat_brain',
        'ncat_spleen',
        'ncat_iodine_blood',
        'ncat_iron',
        'ncat_pmma',
        'ncat_aluminum',
        'ncat_titanium',
        'ncat_air',
        'ncat_graphite',
        'ncat_lead',
        'ncat_breast_mammary',
        'ncat_skin',
        'ncat_iodine',
        'ncat_eye_lens',
        'ncat_ovary',
        'ncat_red_marrow',
        'ncat_yellow_marrow',
        'ncat_testis',
        'ncat_thyroid',
        'ncat_bladder']
    
    Evec = cfg.spec.Evec
    nMat = len(materialList)
    Mus = np.zeros([Evec.size, nMat], dtype=np.double)
    for i in range(nMat):
        Mus[:, i] = GetMu(materialList[i], Evec)/10 # cm^-1 --> mm^-1

    # in nCAT_main.c: void set_material_info_NCAT(int materialCount, int eBinCount, double *muTable)
    # the C func wants data order: materialIndex -> Ebin, so the dim is [Ebin, materialIndex]
    fun = cfg.clib.set_material_info_NCAT
    fun.argtypes = [c_int, c_int, ndpointer(c_double)]
    fun.restype = None
    fun(nMat, Evec.size, Mus)
    
    return nMat


def set_volume(cfg, nMat):
    ###----------- set tolerance, defauts: t1 = 0.001, t2 = 0.0001, pd = 0.02
    # in nCAT_main.c: void set_tolerance_info_NCAT(double t1, double t2, double pd)
    fun_t = cfg.clib.set_tolerance_info_NCAT
    fun_t.argtypes = [c_double, c_double, c_double]
    fun_t.restype = None
    t1 = 0.001
    t2 = 0.0001
    pd = 0.02
    fun_t(t1, t2, pd)
    
    ###----------- set volume
    # phantom file (*.nrb)
    cfg.phantom.filename = my_path.find("phantom", cfg.phantom.filename, '')

    # in nCAT_main.c: void Parse_Phantom(char *filename, int *materials, float *coord_origin_offset, float scale)
    fun = cfg.clib.Parse_Phantom
    fun.argtypes = [POINTER(c_char), ndpointer(c_int), POINTER(c_float), c_float]
    fun.restype = None
    
    phantom_filename = c_char_p(bytes(cfg.phantom.filename, 'utf-8'))
    material_flags = np.ones((1,nMat),dtype=np.int32)
    offsets = (c_float*3)(*cfg.phantom.centerOffset)
    scale = cfg.phantom.scale
    fun(phantom_filename, material_flags, offsets, scale)


def set_detector(cfg):
    det = cfg.det
    
    Height = [det.height]
    Width = [det.width]
    Pix = [det.nCells]
    Coords = det.cellCoords.astype(np.double)
    Sub = [det.nSamples]
    Sampling = det.sampleCoords.astype(np.double)
    Weight = det.weights.astype(np.double)
    nModuleTypes = det.nModDefs
    maxPix = np.max(det.nCells)
    maxSubDets = np.max(det.nSamples)
    UNUSED = 0
    
    Height = (c_double*1)(*Height)
    Width = (c_double*1)(*Width)
    Pix = (c_int*1)(*Pix)
    Sub = (c_int*1)(*Sub)
        
    # in nCAT_main.c: void set_module_info_NCAT(double *Height, double *Width, int *Pix, double *Coords, int *Sub, double *Sampling, double *Weight, int nModuleTypes, int maxPix, int maxSubDets, int UNUSED)
    # the C func wants data order: row -> col -> pixel_ind or sample_ind
    fun = cfg.clib.set_module_info_NCAT
    fun.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_int), ndpointer(c_double), \
        POINTER(c_int), ndpointer(c_double), ndpointer(c_double), c_int, c_int, c_int, c_int]
    fun.restype = None    
    fun(Height, Width, Pix, Coords, \
        Sub, Sampling, Weight, nModuleTypes, maxPix, maxSubDets, UNUSED)
 
 
def set_source(cfg):
    # in nCAT_main.c: void set_src_info_NCAT(double *sourceWeights, int nSubSources)
    fun = cfg.clib.set_src_info_NCAT
    fun.argtypes = [ndpointer(c_double), c_int]
    fun.restype = None
    
    SourceWeights = cfg.src.weights.astype(np.double)
    NumberOfSubSources = cfg.src.nSamples
    
    fun(SourceWeights, NumberOfSubSources)
 
