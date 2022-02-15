import numpy as np
import math


def CreateHSP(Length, index):
    HS = np.ones(Length)
    Center = int((Length) / 2)
    PI = 3.14159265358979

    if index == 1:
        HS[0] = 0
        HS[Center] = 0.25
        for i in range(1, Center):
            HS[i] = -np.power((math.sin(PI * (i - Center) / 2)), 2) / (PI * PI * (i - Center) * (i - Center))

        for k in range((Center + 1), Length):
            HS[k] = -np.power((math.sin(PI * (k - Center) / 2)), 2) / (PI * PI * (k - Center) * (k - Center))

    elif index == 2:
        for i in range(Length):
            HS[i] = -2 / (PI * PI * (4 * (i - Center) * (i - Center) - 1))
            
    else: 
        print("Wrong Kernel type.")

    return HS