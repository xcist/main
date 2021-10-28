import numpy as np
import math

def CreateHSP(Length, index):
    x = 0

    Window = {
        1 : lambda x: np.ones(Length),
        2 : lambda x: np.kaiser(Length, 2.5),
        3 : lambda x: np.hamming(Length),
        4 : lambda x: np.hanning(Length),
        5 : lambda x: np.blackman(51)
    }[index](x)

    HS=np.ones(Length)
    Center = int((Length)/2)
    PI=3.14159265358979
    HS[0]=0
    HS[Center] = 0

    for i in range(1,Center):
        HS[i]=2*np.power((math.sin(PI*(i-Center)/2)),2)/(PI*(i-Center))

    for k in range((Center+1),Length):
        HS[k] = 2*np.power((math.sin(PI*(k-Center)/2)),2)/(PI*(k-Center))

    return HS*Window

