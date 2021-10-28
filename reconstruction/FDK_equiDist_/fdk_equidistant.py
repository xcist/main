import ctypes as ct
import numpy as np
import scipy.io as scio
import matplotlib.pyplot as plt

# Init ctypes types
DOUBLE = ct.c_double
PtrDOUBLE = ct.POINTER(DOUBLE)
PtrPtrDOUBLE = ct.POINTER(PtrDOUBLE)
PtrPtrPtrDOUBLE = ct.POINTER(PtrPtrDOUBLE)


class TestStruct(ct.Structure):
    _fields_ = [
                ("ScanR", ct.c_double),     # Radius of the scanning trajectory
                ("DecWidth", ct.c_double),  # Length of the detector
                ("DecHeigh", ct.c_double),  # z axis length of the detector
                ("YL", ct.c_int),           # Detector cell number on each row
                ("ZL", ct.c_int),           # Detector cell number on each column
                ("YOffSet", ct.c_double),  # Detector offset along the horizontal direction (pixel, e.g. quarter pixel)
                ("ZOffSet", ct.c_double),  # Detector offset along the vertical direcion (pixel, e.g. quarter pixel)
                ("ProjScale", ct.c_int),    # Number of sampling on the scanning trajectory
                ("DistD", ct.c_double),     # Distance between the x-ray source to the detector
                ("Radius", ct.c_double),    # Radius of the phantom
                ("RecSize", ct.c_int),      # Reconstructed image size
                ("centerX", ct.c_int),      # Reconstructed center on x axis
                ("centerY", ct.c_int),      # Reconstructed center on y axis
                ("centerZ", ct.c_int),      # Reconstructed center on z axis
                ("FOILength", ct.c_int),    # Reconstructed length on x axis(pixel number)
                ("FOIWidth", ct.c_int),     # Reconstructed length on y axis(pixel number)
                ("FOIHeigh", ct.c_int),     # Reconstructed length on z axis(pixel number)
                ("GF", PtrPtrPtrDOUBLE),    # Projection data/ Sinogram data
                ("RecIm", PtrPtrPtrDOUBLE)  # Reconstruction data
                ]


def double3darray2pointer(arr):
    # Converts a 3D numpy to ctypes 3D array.
    arr_dimx = DOUBLE * arr.shape[2]
    arr_dimy = PtrDOUBLE * arr.shape[1]
    arr_dimz = PtrPtrDOUBLE * arr.shape[0]

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


# Load the compiled library
recon = ct.CDLL("./fdk_equidist.dll")
# Define arguments of the C function
recon.fbp.argtypes = [ct.POINTER(TestStruct)]
# Define the return type of the C function
recon.fbp.restype = None


# Load the data
dataFile = './data/FDK_Filtering_Dist.mat'
data = scio.loadmat(dataFile)

t = TestStruct()

t.ScanR = data['ScanR']
t.DecWidth = data['DecWidth']
t.DecHeigh = data['DecHeigh']
t.YL = data['YL']
t.ZL = data['ZL']
t.ProjScale = data['ProjScale']
t.DistD = data['DistD']
t.Radius = data['Radius']


# These are flexible parameters.
t.RecSize = 128
t.centerX = 64
t.centerY = 64
t.centerZ = 64
t.FOILength = 128
t.FOIWidth = 128
t.FOIHeigh = 128

# Generate a 2D ctypes array from numpy array
GF = data['GF']
GF_ptr = double3darray2pointer(GF)
t.GF = GF_ptr

RecIm = np.zeros(shape=(t.FOILength, t.FOIWidth, t.FOIHeigh))
RecIm_ptr = double3darray2pointer(RecIm)
t.RecIm = RecIm_ptr

recon.fbp(ct.byref(t))

RecA = double3dpointer2array(RecIm_ptr, *RecIm.shape)

dataNew = './data/FDK_RecImage.mat'
scio.savemat(dataNew,
             {
                 'Rec': RecA,
             },)

plt.figure()
plt.imshow(RecA[:, :, 49], cmap='gray')
plt.show()