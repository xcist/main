import numpy as np
import math
PI = 3.14159265358979


def parker(data, projection_range):
    ScanR = data['ScanR']
    Radius = data['ObjR']
    ProjData = data['Proj']
    ProjScale = int(data['ProjScale'])
    DecFanAng = data['DecAngle']
    Dgy = np.array(ProjData, dtype=np.float32)
    YL = int(data['YL'])

    delta = 2*math.asin(Radius/ScanR)
    dgamma = DecFanAng/YL
    detcenter = (YL+ 1)*0.5
    dbeta = 2*PI*projection_range/(360*ProjScale)

    for pindex in range(ProjScale):
        beta = pindex*dbeta
        for yindex in range(YL):
            gamma = (yindex - detcenter)*dgamma
            if (beta < (delta - 2 * gamma)):
                Dgy[yindex,:,pindex] = (np.sin(PI*beta/(2*delta-4*gamma)))**2
            elif(beta < (PI - 2 * gamma)):
                Dgy[yindex, :, pindex] = 1
            elif (beta < (PI + delta)):
                Dgy[yindex,:,pindex] = (np.sin(PI*(PI + delta - beta)/(2*delta + 4*gamma)))**2
            else:
                Dgy[yindex, :, pindex] = 0
    proj = np.multiply(ProjData, Dgy)

    return proj, Dgy