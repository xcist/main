import ctypes as ct
import numpy as np
import scipy.io as scio
import matplotlib.pyplot as plt


# Init ctypes types
DOUBLE = ct.c_double
PtrDOUBLE = ct.POINTER(DOUBLE)
PtrPtrDOUBLE = ct.POINTER(PtrDOUBLE)


class TestStruct(ct.Structure):
    _fields_ = [
                ("ScanR", ct.c_double),
                ("DecLength", ct.c_double),
                ("YL", ct.c_int),
                ("AngleNumber", ct.c_int),
                ("DistD", ct.c_double),
                ("Radius", ct.c_double),
                ("RecSize", ct.c_int),
                ("centerX", ct.c_int),
                ("centerY", ct.c_int),
                ("FOILength", ct.c_int),
                ("FOIWidth", ct.c_int),
                ("GF", PtrPtrDOUBLE),
                ("RecIm", PtrPtrDOUBLE)
                ]


def double2darray2pointer(array):
    # Converts a 2D numpy into a array ctypes 2D array.
    arr_dimx = DOUBLE * array.shape[1]
    arr_dimy = PtrDOUBLE * array.shape[0]
    arr_ptr = arr_dimy()

    for i, row in enumerate(array):
        arr_ptr[i] = arr_dimx()
        for j, val in enumerate(row):
            arr_ptr[i][j] = val

    return arr_ptr


def double2dpointer2array(ptr, n, m):
    # Converts ctypes 2D array into a 2D numpy array.
    arr = np.zeros(shape=(n, m))
    for i in range(n):
        for j in range(m):
            arr[i][j] = ptr[i][j]
    return arr


# Load the compiled library
recon = ct.CDLL("./fbpequidist.dll")
# Define arguments of the C function
recon.fbp.argtypes = [ct.POINTER(TestStruct)]
# Define the return type of the C function
recon.fbp.restype = None

# Load the data
dataFile = './data/Res_Filtering_Dist.mat'
data = scio.loadmat(dataFile)

# init the struct
t = TestStruct()

t.ScanR = data['ScanR']
t.DecLength = data['DecLength']
t.YL = data['YL']
t.AngleNumber = len(data['GF'])
t.DistD = data['DistD']
t.Radius = data['Radius']

# These are flexible parameters.
t.RecSize = 256
t.centerX = 128
t.centerY = 128
t.FOILength = 256
t.FOIWidth = 256

# Generate a 2D ctypes array from numpy array
GF = data['GF']
GF = GF.T
GF_ptr = double2darray2pointer(GF)
t.GF = GF_ptr

RecIm = np.zeros(shape=(t.FOILength, t.FOIWidth))
RecIm_ptr = double2darray2pointer(RecIm)
t.RecIm = RecIm_ptr

# interface with C function
recon.fbp(ct.byref(t))

# Convert ctypes 2D arrays to numpy arrays
RecA = double2dpointer2array(RecIm_ptr, *RecIm.shape)

dataNew = './data/fbp_equidistant_RecImage.mat'
scio.savemat(dataNew,
             {
                 'Rec': RecA,
             },)

plt.figure()
plt.imshow(RecA, cmap='gray', vmin=0.95, vmax=1.05)
plt.show()
