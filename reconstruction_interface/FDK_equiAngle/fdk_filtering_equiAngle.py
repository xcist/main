import numpy as np
import math
import scipy.io as scio
from CreateHSP import CreateHSP


dataFile = './data/FDK_proj_curve.mat'
data = scio.loadmat(dataFile)

ScanR = data['ScanR']
DistD = data['StdDis']
Radius = data['ObjR']
ProjData = data['Proj']
ProjScale = int(data['ProjScale'])
DecFanAng = data['DecAngle']
Dgy = np.array(ProjData, dtype=np.float32)
YL = int(data['YL'])
ZL = int(data['ZL'])
# Try to import the offset, otherwise set them as zeros
if data.get('YOffSet'):
    YOffSet = data['YOffSet']
else:
    YOffSet = 0

if data.get('ZOffSet'):
    ZOffSet = data['ZOffSet']
else:
    ZOffSet = 0


DecHeigh = data['DecHeigh']
DeltaUW = DecFanAng/(YL-1)
DeltaU2 = 2*DeltaUW

# pre-weighting

for Yindex in range(1, YL-1):
   Dgy[Yindex,:,:]=(ProjData[Yindex+1,:,:]-ProjData[Yindex-1,:,:])/DeltaU2

Dgy[0,:,:] = Dgy[1,:,:]
Dgy[YL-1,:,:]= Dgy[YL-2,:,:]

Dg=Dgy

# filtering

WindowType=1
nn=int(math.pow(2,(math.ceil(math.log2(abs(YL)))+1)))
HS=CreateHSP(nn,WindowType)
nn2= nn*2
k = int(nn/2)
TempF=np.zeros(nn2)
TempF[0:k]=HS[k:nn]
TempF[k+nn:nn2]=HS[0:k]
HS=TempF*complex(0,1)
FFT_F=np.fft.fft(HS)

GF=Dg

for ProjIndex in range(0,ProjScale):
    for j in range(ZL):
       TempData=np.ones(YL)
       for k in range(YL):
         TempData[k]=Dg[k,j,ProjIndex]
       FFT_S=np.fft.fft(TempData,nn2)
       TempData=np.fft.ifft(FFT_S*FFT_F).imag
       for k in range(YL):
        GF[k,j,ProjIndex]=-TempData[k]


dataNew = './data/FDK_Filtering_curve.mat'
scio.savemat(dataNew,
             {'GF': Dgy,
             'ScanR': ScanR,
             'DistD': DistD,
             'DecFanAng': DecFanAng,
             'ProjScale': ProjScale,
             'YL': YL,
             'YOffSet': YOffSet,
             'DecHeigh': DecHeigh,
             'ZL': ZL,
             'ZOffSet':ZOffSet,
             'Radius': Radius,},
             )
