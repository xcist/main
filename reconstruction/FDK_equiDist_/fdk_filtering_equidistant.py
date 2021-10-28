import numpy as np
import math
import scipy.io as scio
from CreateHSP import CreateHSP


# dataFile = './data/FDK_proj.mat'
dataFile = './data/ProjRes.mat'
data = scio.loadmat(dataFile)

ScanR = data['ScanR']
DistD = int(data['StdDis'])
Radius = data['ObjR']
ProjScale = int(data['ProjScale'])
ProjData = data['Proj']
Dgu = np.array(ProjData, dtype=np.float32)
YL = int(data['YL'])
ZL = int(data['ZL'])
HalfY = (YL+1)*0.5
HalfZ = (ZL+1)*0.5
DecWidth = data['DecWidth']
DecHeigh = data['DecHeigh']
DeltaUW = DecWidth/data['YL']
DeltaU2 = 2*DeltaUW
DecWidth = data['DecWidth']
DeltaUH = DecHeigh/data['ZL']


if data.get('YOffSet'):
    YOffSet = data['YOffSet']
else:
    YOffSet = 0

if data.get('ZOffSet'):
    ZOffSet = data['ZOffSet']
else:
    ZOffSet = 0

# pre-weighting

W = np.zeros((YL,ZL))
for Zindex in range(1, ZL-1):
    for Yindex in range(1, YL-1):
        W[Yindex,Zindex] = np.sqrt(DistD * DistD + np.power((Yindex - HalfY) * DeltaUW, 2)+ np.power((Zindex - HalfZ) * DeltaUH, 2)) / DistD
        Dgu[Yindex,Zindex,:]=W[Yindex, Zindex]*(ProjData[Yindex+1,Zindex,:]-ProjData[Yindex-1,Zindex,:])/DeltaU2

Dgu[0,:,:] = Dgu[1,:,:]
Dgu[YL-1,:,:] = Dgu[YL-2,:,:]
Dgu[:,0,:]= Dgu[:,1,:]
Dgu[:,ZL-1,:]= Dgu[:,ZL-2,:]

Dg=Dgu

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
        GF[k,j,ProjIndex]=TempData[k]


dataNew = './data/FDK_Filtering_Dist.mat'
scio.savemat(dataNew,
             {'GF': GF,
             'ScanR': ScanR,
             'DistD': DistD,
             'ProjScale': ProjScale,
             'DecWidth': DecWidth,
             'YL': YL,
             'YOffSet': YOffSet,
             'DecHeigh': DecHeigh,
             'ZL': ZL,
             'ZOffSet': ZOffSet,
             'Radius': Radius,},
             )
