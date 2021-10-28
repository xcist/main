import numpy as np
import math
import scipy.io as scio
from CreateHSP import CreateHSP
import matplotlib.pyplot as plt


dataFile = './data/Res_Proj.mat'
data = scio.loadmat(dataFile)

Number = len(data['ProjData'])
ProjData = data['ProjData'][0:Number,:]
Angle = data['Angle'][0,:]
YL = int(data['YL'])
HalfY = (YL+1)*0.5
DecLength = data['DecLength']
DeltaU = DecLength/data['YL']
DeltaU2 = 2*DeltaU
Dgu = np.array(ProjData, dtype=np.float32)
DistD = int(data['DistD'])
Radius = data['Radius']
ScanR = data['ScanR']
W = np.zeros(YL)

for Yindex in range(1, YL-1):
    W[Yindex] = np.sqrt(DistD * DistD + np.power((Yindex - HalfY) * DeltaU, 2)) / DistD
    Dgu[:,Yindex]=W[Yindex]*(ProjData[:,Yindex+1]-ProjData[:,Yindex-1])/DeltaU2

Dgu[:,0] = Dgu[:,1]
Dgu[:,YL-1]= Dgu[:,YL-2]
Dg=Dgu


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
for ProjIndex in range(0,Number):
   TempData=np.ones(YL)
   for k in range(0, YL):
     TempData[k]=Dg[ProjIndex,k]
   FFT_S=np.fft.fft(TempData,nn2)
   TempData=np.fft.ifft(FFT_S*FFT_F).imag
   for k in range(0,YL):
    GF[ProjIndex,k]=TempData[k]

GF = GF[1:Number-1,:]
Angle = Angle[1:Number-1]

dataNew = './data/Res_Filtering_Dist.mat'
scio.savemat(dataNew,
             {'GF': GF,
             'Angle': Angle,
             'Radius': Radius,
             'ScanR': ScanR,
             'DistD': DistD,
             'DecLength': DecLength,
             'YL': YL},
             )
