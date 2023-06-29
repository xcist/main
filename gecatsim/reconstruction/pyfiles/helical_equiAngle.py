import ctypes as ct
import numpy as np
import math
from gecatsim.reconstruction.pyfiles.createHSP import createHSP
import scipy.io as scio
import matplotlib.pyplot as plt
from gecatsim.reconstruction.pyfiles.mapConfigVariablesToHelical import mapConfigVariablesToHelical
import os

# Init ctypes types
FLOAT = ct.c_float
PtrFLOAT = ct.POINTER(FLOAT)
PtrPtrFLOAT = ct.POINTER(PtrFLOAT)
PtrPtrPtrFLOAT = ct.POINTER(PtrPtrFLOAT)
PI =3.14159265358979

class TestStruct(ct.Structure):
    _fields_ = [
                ("ScanR", ct.c_float),
                ("DistD", ct.c_float),
                ("YL", ct.c_int),
                ("ZL", ct.c_int),
                ("dectorYoffset", ct.c_float),  # Detector along the horizontal direction (pixel, e.g. quarter pixel)
                ("dectorZoffset", ct.c_float),  # Detector offset along the vertical direcion (pixel, e.g. quarter pixel)
                ("XOffSet", ct.c_float),  # recon offset along the x axis(mm)
                ("YOffSet", ct.c_float),  # recon offset along the y axis(mm)
                ("ZOffSet", ct.c_float),  # recon offset along the z axis(mm)
                ("phantomXOffSet", ct.c_float),  # phantom offset along the x axis(mm)
                ("phantomYOffSet", ct.c_float),  # phantom offset along the y axis(mm)
                ("phantomZOffSet", ct.c_float),  # phantom offset along the z axis(mm)
                ("DecFanAng", ct.c_float),
                ("DecHeight", ct.c_float),
                ("DecWidth", ct.c_float),
                ("dx", ct.c_float),
                ("dy", ct.c_float),
                ("dz", ct.c_float),
                ("h", ct.c_float),
                ("BetaS", ct.c_float),
                ("BetaE", ct.c_float),
                ("AngleNumber", ct.c_int),
                ("N_2pi", ct.c_int),
                ("Radius", ct.c_float),
                ("RecSize", ct.c_int),
                ("RecSizeZ", ct.c_int),
                ("delta", ct.c_float),
                ("HSCoef", ct.c_float),
                ("k1", ct.c_float),
                ("GF", PtrPtrPtrFLOAT),
                ("RecIm", PtrPtrPtrFLOAT)
                ]


def double3darray2pointer(arr):
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


def double3dpointer2array(ptr, n, m, o):
    # Converts ctypes 3D array into a 3D numpy array.
    arr = np.zeros(shape=(n, m, o))

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
        lib_file = "helicalrecon.dll"
    else:
        lib_file = "helicalrecon.so"
    clib = ll(os.path.join(recon_lib, lib_file))

    return clib


def helical_equiAngle(cfg, prep):
    prep = prep[:,:,::-1]

    if cfg.protocol.rotationDirection == -1:	
        prep = prep[::-1,:,:]
    Proj = prep.transpose(0,2,1)

    # scanner & recon geometry
    SO, DD, YL, ZL, ViewN, N_Turn, N_2pi, h, startAngle, dectorYoffset, dectorZoffset, \
    k1, delta, HSCoef, rowSize, ColSize, imageSize, sliceCount, sliceThickness, centerOffset, ObjR,kernelType \
        = mapConfigVariablesToHelical(cfg)

    # print(SO, DD, YL, ZL, ViewN, N_Turn, N_2pi, h, dectorYoffset, dectorZoffset, imageSize, sliceCount, ObjR, ColSize)
    
    DecAngle = math.atan(ColSize/2/DD)*2*YL
    HSCoef = (DecAngle/PI+0.52)
    DecHeight = rowSize * ZL
    h = h/DecHeight
    dx = 2*ObjR/imageSize
    dy = 2*ObjR/imageSize
    dz = sliceThickness

    YLCW = (YL-1)*0.5
    YLC= (YL-1)*0.5 - dectorYoffset  # Detector center along the horizontal direction of detector array
    ZLC= (ZL-1)*0.5 - dectorZoffset # Detector center along the vertical direction of detector array
    
    BetaS = -N_Turn*PI  + PI*startAngle/180
    BetaE =  N_Turn*PI  + PI*startAngle/180

    DecWidth = math.tan(DecAngle*0.5)*(SO)*2
    dYL   =  DecWidth/YL
    dZL   =  DecHeight/ZL
    DeltaFai= (BetaE-BetaS)/(ViewN-1)
    DeltaTheta = DeltaFai
    DeltaT     = dYL  
    dYA = DecAngle/YL
    PProj    = np.zeros((ViewN,YL,ZL))

    ## rebinning the projection
    print("* Rebinning the projection...")
    for i in range(ViewN):
        Theta=(i)*DeltaTheta                   # the view for the parallel projection
        for j in range(YL):
            t      = (j-YLCW)*DeltaT     # the distance from origin to ray for parallel beam
            Beta   = math.asin(t/(SO))            # the fan_angle for cone_beam projection
            Fai    = Theta+Beta              # the view for cone_beam projecton
            a      = Beta  # the position of this ray on the flat detector
            FaiIndex        =  (Fai/DeltaFai)
            UIndex          =  (a/dYA)+YLC
            FI              =  math.ceil(FaiIndex)
            UI              =  math.ceil(UIndex)
            coeXB           =  FI-FaiIndex
            coeXU           =  1-coeXB
            coeYB           =  UI-UIndex
            coeYU           =  1-coeYB
            if (FI<=0):
                IndexXU = 0
                IndexXB = 0
            elif(FI > ViewN-1):
                IndexXU = ViewN-1
                IndexXB = ViewN-1
            else:
                IndexXU = FI
                IndexXB = FI-1

            if (UI<=0):
                IndexYU = 0
                IndexYB = 0
            elif(UI>YL-1):
                IndexYU = YL-1
                IndexYB = YL-1
            else:
                IndexYU=UI
                IndexYB=UI-1
            PProj[i,j,:]=coeXB*coeYB*Proj[IndexXB,IndexYB,:]+coeXU*coeYB*Proj[IndexXU,IndexYB,:]+coeXB*coeYU*Proj[IndexXB,IndexYU,:]+coeXU*coeYU*Proj[IndexXU,IndexYU,:]

    Proj = PProj.transpose(1,2,0)

    #scio.savemat('testrebin.mat', {'rebin': PProj})

    #Perform Ramp filtering
    print("* Applying the filter...")
    Dg=Proj
    nn = int(math.pow(2, (math.ceil(math.log2(abs(YL))) + 1)))
    nn2 = nn*2
    FFT_F = createHSP(nn, kernelType)

    GF = Proj
    
    for ProjIndex in range(0, ViewN):
        for j in range(ZL):
            TempData = np.ones(YL)
            for k in range(YL):
                TempData[k] = Dg[k, j, ProjIndex]
            FFT_S = np.fft.fft(TempData, nn2)
            TempData = np.fft.ifft(FFT_S * FFT_F).imag
            for k in range(YL):
                GF[k, j, ProjIndex] = TempData[k]
    GF = GF/dYL

    #Backproject the filtered data into the 3D space
    # Load the compiled library
    recon = load_C_recon_lib()
    # Define arguments of the C function
    recon.fbp.argtypes = [ct.POINTER(TestStruct)]
    # Define the return type of the C function
    recon.fbp.restype = None

    # init the struct
    t = TestStruct()

    t.ScanR = SO
    t.DistD = DD
    t.YL = YL
    t.ZL = ZL
    t.DecFanAng = DecAngle
    t.DecHeight = DecHeight
    t.DecWidth = DecWidth
    t.dx = dx
    t.dy = dy 
    t.dz = dz
    t.h = h
    t.BetaS = BetaS
    t.BetaE = BetaE
    t.dectorYoffset = dectorYoffset
    t.dectorZoffset = dectorZoffset
    t.AngleNumber = ViewN
    t.N_2pi = N_2pi
    t.Radius = ObjR
    t.RecSize = imageSize
    t.RecSizeZ = sliceCount
    t.delta = delta
    t.HSCoef = HSCoef
    t.k1 = k1

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
        # print("* Start view: {}".format(t.startAngle))
        print("* Number of detector cols: {}".format(t.YL))
        print("* Number of detector rows: {}".format(t.ZL))
        print("* Detector height: {} mm".format(t.DecHeight))
        print("* Detector X offset: {} mm".format(t.dectorYoffset))
        print("* Detector Z offset: {} mm".format(t.dectorZoffset))
        print("* Scan number of views: {} ".format(t.AngleNumber))
        print("* Recon FOV: {} mm".format(2 * t.Radius))
        print("* Recon XY pixel size: {} mm".format(t.RecSize))
        print("* Recon Slice thickness: {} mm".format(t.sliceThickness))
        print("* Recon X offset: {} mm".format(t.XOffSet))
        print("* Recon Y offset: {} mm".format(t.YOffSet))
        print("* Recon Z offset: {} mm".format(t.ZOffSet))
    # Generate a 3D ctypes array from numpy array
    print("* Converting projection data from a numpy array to a C array...")
    GF_ptr = double3darray2pointer(GF)
    t.GF = GF_ptr

    # RecIm = np.zeros(shape=(t.RecSize, t.RecSize, t.RecSize))
    print("* Allocating a C array for the recon results...")
    RecIm = np.zeros(shape=(t.RecSize, t.RecSize, t.RecSizeZ))
    RecIm_ptr = double3darray2pointer(RecIm)
    t.RecIm = RecIm_ptr

    # interface with C function
    print("* In C...")
    recon.fbp(ct.byref(t))

    # Convert ctypes 3D arrays to numpy arrays
    print("* Converting the recon results from a C array to a numpy array...")
    rec = double3dpointer2array(RecIm_ptr, *RecIm.shape)
    rec = rec.transpose(1,0,2).astype(np.float32)
    rec = rec[::-1,::-1]

    return rec
