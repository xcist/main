import numpy as np
import math
import scipy.io as scio
from CreateHSP import CreateHSP
import matplotlib.pyplot as plt

dataFile = './data/Res_Proj_Angle.mat'
data = scio.loadmat(dataFile)

Number = len(data['ProjData'])
ProjData = data['ProjData'][0:Number,:]
Angle = data['Angle'][0,:]
YL = int(data['YL'])
DecFanAng = data['DecFanAng']
DeltaY = DecFanAng/(YL-1)
DeltaY2 = 2*DeltaY
Dgy = np.array(ProjData, dtype=np.float32)
Radius = data['Radius']
ScanR = data['ScanR']

for Yindex in range(1, YL-1):
   Dgy[:,Yindex]=(ProjData[:,Yindex+1]-ProjData[:,Yindex-1])/DeltaY

Dgy[:,0] = Dgy[:,1]
Dgy[:,YL-1]= Dgy[:,YL-2]


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

GF=Dgy
for ProjIndex in range(0,Number):
   TempData=np.ones(YL)
   for k in range(0, YL):
     TempData[k]=Dgy[ProjIndex,k]
   FFT_S=np.fft.fft(TempData,nn2)
   TempData=np.fft.ifft(FFT_S*FFT_F).imag
   for k in range(0,YL):
    GF[ProjIndex,k]=-TempData[k]

GF = GF[0:Number,:]
Angle = Angle[0:Number]


dataNew = './data/Res_Filtering_Angle.mat'
scio.savemat(dataNew,
             {'GF': GF,
             'Angle': Angle,
             'Radius': Radius,
             'ScanR': ScanR,
             'DecFanAng': DecFanAng,
             'YL': YL},
             )