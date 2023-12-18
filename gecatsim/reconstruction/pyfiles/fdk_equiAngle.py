# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import ctypes as ct
import math
from gecatsim.reconstruction.pyfiles.createHSP import createHSP
from gecatsim.reconstruction.pyfiles.mapConfigVariablesToFDK import mapConfigVariablesToFDK
from gecatsim.pyfiles.CommonTools import *

# Init ctypes types
FLOAT = ct.c_float
PtrFLOAT = ct.POINTER(FLOAT)
PtrPtrFLOAT = ct.POINTER(PtrFLOAT)
PtrPtrPtrFLOAT = ct.POINTER(PtrPtrFLOAT)


class TestStruct(ct.Structure):
    _fields_ = [
                ("ScanR", ct.c_float),         # Radius of the scanning trajectory of x-ray source
                ("DecFanAng", ct.c_float),     # Fan angle coverage of the detector element along the horizontal diretion
                ("startangle", ct.c_float),     # recon startangle
                ("rotdir", ct.c_int), # //rotation direction
                ("DecHeight", ct.c_float),      # Physical height of the detector along the vertical direction
                ("YL", ct.c_int),               # Detector cell number on each row along the horizontal direction
                ("ZL", ct.c_int),               # Detector cell number on each column along the vertical direction
                ("dectorYoffset", ct.c_float),               # Detector along the horizontal direction (pixel, e.g. quarter pixel)
                ("dectorZoffset", ct.c_float),               # Detector offset along the vertical direcion (pixel, e.g. quarter pixel)
                ("XOffSet", ct.c_float),       # recon offset along the x axis
                ("YOffSet", ct.c_float),       # recon offset along the y axis
                ("ZOffSet", ct.c_float),       # recon offset along the z axis
                ("phantomXOffSet", ct.c_float),  # phantom offset along the x axis
                ("phantomYOffSet", ct.c_float),  # phantom offset along the y axis
                ("phantomZOffSet", ct.c_float),  # phantom offset along the z axis
                ("AngleNumber", ct.c_int),      # Number of view samples on the scanning trajectory
                ("DistD", ct.c_float),         # Distance between the x-ray source and the detector
                ("Radius", ct.c_float),        # Radius of the phantom
                ("RecSize", ct.c_int),          # Reconstructed size
                ("sliceThickness", ct.c_float),   # Reconstructed thickness
#                 ("centerX", ct.c_float),          # Reconstructed center on x axis
#                 ("centerY", ct.c_float),          # Reconstructed center on y axis
#                 ("centerZ", ct.c_float),          # Reconstructed center on z axis
                ("FOILength", ct.c_int),        # Reconstructed length on x axis
                ("FOIWidth", ct.c_int),         # Reconstructed length on y axis
                ("FOIHeight", ct.c_int),         # Reconstructed length on z axis
                ("GF", PtrPtrPtrFLOAT),        # Projection data/ Sinogram data
                ("RecIm", PtrPtrPtrFLOAT)      # Reconstructed 3D image data
                ]


def float3Darray2pointer(arr):
    # Converts a 3D numpy to ctypes 3D array.
    arr_dimx = FLOAT * arr.shape[2]
    arr_dimy = PtrFLOAT * arr.shape[1]
    arr_dimz = PtrPtrFLOAT * arr.shape[0]

    arr_ptr = arr_dimz()

    for i, row in enumerate(arr):
        arr_ptr[i] = arr_dimy()
        for j, col in enumerate(row):
            arr_ptr[i][j] = arr_dimx()
            for k, val in enumerate(col):
                arr_ptr[i][j][k] = val
    return arr_ptr


def float3Dpointer2array(ptr, n, m, o):
    # Converts ctypes 3D array into a 3D numpy array.
    arr = np.zeros(shape=(n, m, o), dtype=np.single)

    for i in range(n):
        for j in range(m):
            for k in range(o):
                arr[i, j, k] = ptr[i][j][k]

    return arr


def load_C_recon_lib():

    # add recon lib path to environment value "PATH" for depending DLLs
    # # # # recon_lib = my_path.find_dir("top", os.path.join("reconstruction", "lib"))
    # # # # my_path.add_dir_to_path(recon_lib)
    
    #  my_path.find_dir doesn't have the key "reconstruction", use the temp solution below:
    recon_lib = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib")
    

    # load C/C++ lib
    ll = ct.cdll.LoadLibrary
    if os.name == "nt":
        lib_file = "fdk_equiAngle.dll"
    elif os.uname().sysname == 'Darwin':
        lib_file = "fdk_equiAngle_macos.so"
    else:
        lib_file = "fdk_equiAngle.so"
    clib = ll(os.path.join(recon_lib, lib_file))
    
    return clib


def fdk_equiAngle(cfg, prep):

    # scanner & recon geometry
    sid, sdd, nMod, rowSize, modWidth, dectorYoffset, dectorZoffset, \
    fov, imageSize, sliceCount, sliceThickness, centerOffset, startView, rotdir, kernelType      \
        = mapConfigVariablesToFDK(cfg)

    # initialize parameters
    ScanR = sid
    DistD = sdd
    Radius = fov/2
    ProjData = prep.transpose(2,1,0)
    # ProjData = ProjData[::-1,:,:]
    ProjScale = cfg.protocol.viewsPerRotation
    DecFanAng = nMod*2*math.atan(modWidth/2/sdd)
    Dgy = np.array(ProjData, dtype=np.float32)
    YL = int(cfg.scanner.detectorColCount)
    ZL = int(cfg.scanner.detectorRowCount)

    YCtr = (YL - 1) * 0.5
    ZCtr = (ZL - 1) * 0.5
    DecHeight = rowSize*ZL
    
    DeltaUW = DecFanAng/ YL
    DeltaZ = DecHeight / ZL
    
    
    ############## pre-weighting for ramp-filter

    print("* Pre-weighting the filter...")
    for Yindex in range(YL):
        for zindex in range(ZL):
            Dgy[Yindex, zindex, :] = (DistD / np.sqrt(DistD ** 2 + ((zindex - ZCtr) * DeltaZ) ** 2)) * ProjData[Yindex,
                                    zindex,:] * math.cos((Yindex - YCtr) * DeltaUW)

    Dg=Dgy

    ############## filtering

    print("* Applying the filter...")
    nn = int(math.pow(2, (math.ceil(math.log2(abs(YL))) + 1)))
    nn2 = nn*2
    FFT_F = createHSP(nn, kernelType)

    GF = Dg

    for ProjIndex in range(0, ProjScale):
        for j in range(ZL):
            TempData = np.ones(YL)
            for k in range(YL):
                TempData[k] = Dg[k, j, ProjIndex]
            FFT_S = np.fft.fft(TempData, nn2)
            TempData = np.fft.ifft(FFT_S * FFT_F).imag
            for k in range(YL):
                GF[k, j, ProjIndex] = -TempData[k]
    
    GF = GF/DeltaUW
    
    # special case when ZL is 1
    if ZL == 1:
        GF = np.append(GF, GF, axis=1)
        ZL = 2


    ############## FBP
 
    print("* Running the reconstruction...")
    # Load the compiled library
    recon = load_C_recon_lib()

    # Define arguments of the C function
    recon.fbp.argtypes = [ct.POINTER(TestStruct)]

    # Define the return type of the C function
    recon.fbp.restype = None
    
    # init the struct
    t = TestStruct()

    t.ScanR = ScanR
    t.DistD = DistD
    t.DecFanAng = DecFanAng
    t.startangle = startView # Inconsistent C variable name - should be startAngle. Also, this seems to be the start view, not angle.
    t.rotdir = rotdir
    t.DecHeight = DecHeight
    t.YL = YL
    t.ZL = ZL
    t.dectorYoffset = dectorYoffset
    t.dectorZoffset = dectorZoffset
    t.AngleNumber = ProjScale

    t.Radius = Radius
    t.RecSize = imageSize
    t.sliceThickness = sliceThickness
    t.FOILength = imageSize
    t.FOIWidth = imageSize
    t.FOIHeight = sliceCount
#     t.centerX = (t.RecSize - 1)*0.5
#     t.centerY = (t.RecSize - 1)*0.5
#     t.centerZ = (t.FOIHeight - 1)*0.5
    t.XOffSet = centerOffset[0]
    t.YOffSet = centerOffset[1]
    t.ZOffSet = centerOffset[2]
    t.phantomXOffSet = 0
    t.phantomYOffSet = 0
    t.phantomZOffSet = 0
    
    if cfg.recon.printReconParameters:
        print("* Reconstruction parameters:")
        print("* SID: {} mm".format(t.ScanR))
        print("* SDD: {} mm".format(t.DistD))
        print("* Fan angle: {} degrees".format(t.DecFanAng))
        print("* Start view: {}".format(t.startangle))
        print("* Number of detector cols: {}".format(t.YL))
        print("* Number of detector rows: {}".format(t.ZL))
        print("* Detector height: {} mm".format(t.DecHeight))
        print("* Detector X offset: {} mm".format(t.dectorYoffset))
        print("* Detector Z offset: {} mm".format(t.dectorZoffset))
        print("* Scan number of views: {} ".format(t.AngleNumber))
        print("* Recon FOV: {} mm".format(2*t.Radius))
        print("* Recon XY pixel size: {} mm".format(t.RecSize))
        print("* Recon Slice thickness: {} mm".format(t.sliceThickness))
        print("* Recon XY: {} pixels".format(t.FOIWidth))
        print("* Recon Z: {} slices".format(t.FOIHeight))
        print("* Recon X offset: {} mm".format(t.XOffSet))
        print("* Recon Y offset: {} mm".format(t.YOffSet))
        print("* Recon Z offset: {} mm".format(t.ZOffSet))

    print("* Converting projection data from a numpy array to a C array...")
    GF_ptr = float3Darray2pointer(GF)
    t.GF = GF_ptr

    print("* Allocating a C array for the recon results...")
    RecIm = np.zeros(shape=(t.FOILength, t.FOIWidth, t.FOIHeight), dtype=np.single)
    RecIm_ptr = float3Darray2pointer(RecIm)
    t.RecIm = RecIm_ptr

    print("* In C...")
    recon.fbp(ct.byref(t))

    print("* Converting the recon results from a C array to a numpy array...")
    RecA = float3Dpointer2array(RecIm_ptr, *RecIm.shape)
    RecA = RecA[:,:,::-1]

    return RecA
