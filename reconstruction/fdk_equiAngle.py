import ctypes as ct
import numpy as np
import math
import scipy.io as scio
import matplotlib.pyplot as plt
from reconstruction.CreateHSP import CreateHSP
from reconstruction.mapConfigVariablesToFDK import mapConfigVariablesToFDK
from catsim.CommonTools import *
import scipy.interpolate
from scipy.interpolate import interp1d


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
                ("DecHeight", ct.c_float),      # Physical height of the detector along the vertical direction
                ("YL", ct.c_int),               # Detector cell number on each row along the horizontal direction
                ("ZL", ct.c_int),               # Detector cell number on each column along the vertical direction
                ("dectorYoffset", ct.c_float),               # Detector along the horizontal direction (pixel, e.g. quarter pixel)
                ("dectorYoffset", ct.c_float),               # Detector offset along the vertical direcion (pixel, e.g. quarter pixel)
                ("XOffSet", ct.c_float),       # recon offset along the x axis
                ("YOffSet", ct.c_float),       # recon offset along the y axis
                ("ZOffSet", ct.c_float),       # recon offset along the z axis
                ("AngleNumber", ct.c_int),      # Number of view samples on the scanning trajectory
                ("DistD", ct.c_float),         # Distance between the x-ray source and the detector
                ("Radius", ct.c_float),        # Radius of the phantom
                ("RecSize", ct.c_int),          # Reconstructed size
                ("sliceThickness", ct.c_float),   # Reconstructed thickness
                ("centerX", ct.c_float),          # Reconstructed center on x axis
                ("centerY", ct.c_float),          # Reconstructed center on y axis
                ("centerZ", ct.c_float),          # Reconstructed center on z axis
                ("FOILength", ct.c_int),        # Reconstructed length on x axis
                ("FOIWidth", ct.c_int),         # Reconstructed length on y axis
                ("FOIHeight", ct.c_int),         # Reconstructed length on z axis
                ("GF", PtrPtrPtrFLOAT),        # Projection data/ Sinogram data
                ("RecIm", PtrPtrPtrFLOAT)      # Reconstructed 3D image data
                ]



def double3darray2pointer(arr):
    print('* In double3darray2pointer')
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
    print('* In double3dpointer2array')
    # Converts ctypes 3D array into a 3D numpy array.
    arr = np.zeros(shape=(n, m, o))

    for i in range(n):
        for j in range(m):
            for k in range(o):
                arr[i, j, k] = ptr[i][j][k]

    return arr


def fdk_equiAngle(cfg, prep):

    # scanner & recon geometry
    sid, sdd, nRow, nCol, nMod, rowSize, colSize, modWidth, dectorYoffset, dectorZoffset, \
    fov, imageSize, sliceCount, sliceThickness, centerOffset, startAngle, kernelType      \
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
    YOffSet = dectorYoffset
    ZOffSet = 0
    DecHeight = rowSize*ZL
    
    DeltaUW = DecFanAng/(YL-1)
    DeltaU2 = 2*DeltaUW
    DeltaZ = DecHeight / ZL
    
    
    ############## pre-weighting for ramp-filter

    print('* Building the filter')

    for Yindex in range(YL):
        for zindex in range(ZL):
            Dgy[Yindex, zindex, :] = (DistD / np.sqrt(DistD ** 2 + ((zindex - ZCtr) * DeltaZ) ** 2)) * ProjData[Yindex,
                                    zindex,:] * math.cos((Yindex - YCtr) * DeltaUW)

    Dg=Dgy

    ############## filtering

    print('* Applying the filter')

    nn = int(math.pow(2, (math.ceil(math.log2(abs(YL))) + 1)))
    nn2 = nn*2
    FFT_F = CreateHSP(nn, kernelType)

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

    # special case when ZL is 1
    if ZL == 1:
        GF = np.append(GF, GF, axis=1)
        ZL = 2


    ############## FBP
 
    print('* Running the reconstruction *')

    # Load the compiled library
    recon = ct.CDLL("C:/Users/200003237/Documents/GitHub/CatSim/reconstruction/lib/fdk_equiAngle.dll")
    # Define arguments of the C function
    recon.fbp.argtypes = [ct.POINTER(TestStruct)]
    # Define the return type of the C function
    recon.fbp.restype = None
    
    # init the struct
    t = TestStruct()

    t.ScanR = ScanR
    t.DistD = DistD
    t.DecFanAng = DecFanAng
    t.startangle = startAngle + 180
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
    t.centerX = (t.RecSize - 1)*0.5
    t.centerY = (t.RecSize - 1)*0.5
    t.centerZ = (t.FOIHeight - 1)*0.5
    t.XOffSet = centerOffset[0]
    t.YOffSet = centerOffset[1]
    t.ZOffSet = centerOffset[2]

    
    # Generate a 2D ctypes array from numpy array
    GF_ptr = double3darray2pointer(GF)
    t.GF = GF_ptr

    # RecIm = np.zeros(shape=(t.RecSize, t.RecSize, t.RecSize))
    RecIm = np.zeros(shape=(t.FOILength, t.FOIWidth, t.FOIHeight))
    RecIm_ptr = double3darray2pointer(RecIm)
    t.RecIm = RecIm_ptr

    # interface with C function
    print('* In C...')
    recon.fbp(ct.byref(t))

    # Convert ctypes 2D arrays to numpy arrays
    RecA = double3dpointer2array(RecIm_ptr, *RecIm.shape)
    RecA = RecA[:,:,::-1]


    ##--------- Show results

    if t.FOIHeight== 1:
        sliceIndicesToPlot = [0]
    elif t.FOIHeight<= 16:
        sliceIndicesToPlot = range(0, t.FOIHeight)
    else:
        sliceIndicesToPlot = [0, 1, 2,
                              int(t.FOIHeight/2), 
                              t.FOIHeight-3, t.FOIHeight-2, t.FOIHeight-1]

    for sliceIndexToPlot in sliceIndicesToPlot:
        plt.figure(int(sliceIndexToPlot+1))
        sliceToPlot = RecA[:, :, sliceIndexToPlot]
        plt.imshow(sliceToPlot, cmap='gray')
        plt.title("slice " + str(sliceIndexToPlot+1) + " of " + str(t.FOIHeight))
        sliceToPlot = sliceToPlot.copy(order='C')
        fileName = "Slice" + str(sliceIndexToPlot+1) + "_" + str(t.FOILength) + "r" + str(t.FOIWidth) + "s1.raw"
        rawwrite(fileName, sliceToPlot)

    plt.draw()
    plt.pause(1)
    print('*******************************************')
    print('* Press Enter to close plots and continue *')
    input('*******************************************')
    plt.close('all')

    return RecA