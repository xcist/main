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
                ("ScanR", ct.c_double),         # Radius of the scanning trajectory of x-ray source
                ("DecFanAng", ct.c_double),     # Fan angle coverage of the detector element along the horizontal diretion
                ("DecHeigh", ct.c_double),      # Physical heigth of the detector along the vertical direction
                ("YL", ct.c_int),               # Detector cell number on each row along the horizontal direction
                ("ZL", ct.c_int),               # Detector cell number on each column along the vertical direction
                ("YOffSet", ct.c_double),       # Detector offset along the horizontal direction (pixel, e.g. quarter pixel)  
                ("ZOffSet", ct.c_double),       # Detector offset along the vertical direcion (pixel, e.g. quarter pixel)
                ("AngleNumber", ct.c_int),      # Number of view samples on the scanning trajectory
                ("DistD", ct.c_double),         # Distance between the x-ray source and the detector 
                ("Radius", ct.c_double),        # Radius of the phantom
                ("RecSize", ct.c_int),          # Reconstructed size
                ("centerX", ct.c_int),          # Reconstructed center on x axis
                ("centerY", ct.c_int),          # Reconstructed center on y axis
                ("centerZ", ct.c_int),          # Reconstructed center on z axis
                ("FOILength", ct.c_int),        # Reconstructed length on x axis
                ("FOIWidth", ct.c_int),         # Reconstructed length on y axis
                ("FOIHeigh", ct.c_int),         # Reconstructed length on z axis
                ("GF", PtrPtrPtrDOUBLE),        # Projection data/ Sinogram data
                ("RecIm", PtrPtrPtrDOUBLE)      # Reconstructed 3D image data
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
recon = ct.CDLL("./fdk_equiAngle.dll")
# Define arguments of the C function
recon.fbp.argtypes = [ct.POINTER(TestStruct)]
# Define the return type of the C function
recon.fbp.restype = None


# Load the data
dataFile = './data/FDK_Filtering_curve.mat'
data = scio.loadmat(dataFile)

# init the struct
t = TestStruct()

t.ScanR = data['ScanR']
t.DistD = data['DistD']
t.DecFanAng = data['DecFanAng']
t.DecHeigh = data['DecHeigh']
t.YL = data['YL']
t.ZL = data['ZL']

t.AngleNumber = data['ProjScale']
t.Radius = data['Radius']

# These are flexible parameters.
t.RecSize = 128
t.centerX = 64
t.centerY = 64
t.centerZ = 64
t.FOILength = 128
t.FOIWidth = 128
t.FOIHeigh = 1


# Generate a 2D ctypes array from numpy array
GF = data['GF']
GF_ptr = double3darray2pointer(GF)
t.GF = GF_ptr

# RecIm = np.zeros(shape=(t.RecSize, t.RecSize, t.RecSize))
RecIm = np.zeros(shape=(t.FOILength, t.FOIWidth, t.FOIHeigh))
RecIm_ptr = double3darray2pointer(RecIm)
t.RecIm = RecIm_ptr

# interface with C function
recon.fbp(ct.byref(t))

# Convert ctypes 2D arrays to numpy arrays
RecA = double3dpointer2array(RecIm_ptr, *RecIm.shape)

# save result
dataNew = './data/FDK_RecImage_curve.mat'
scio.savemat(dataNew,
             {'Rec': RecA})

plt.figure()
plt.imshow(RecA[:, :, 0], cmap='gray')
plt.show()