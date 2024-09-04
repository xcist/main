# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
from gecatsim.pyfiles.GetMu import GetMu


Mu = []
Mu.append(GetMu('water', 70))
Mu.append(GetMu('water', 70.0))
Mu.append(GetMu('bone', (30, 50, 70)))
Mu.append(GetMu('bone', [30, 50, 70]))
Mu.append(GetMu('Al', range(10, 60, 10)))
types = ['int', 'float', 'tuple', 'list', 'range']

for mu in Mu:
    print()
    for ii in range(len(mu)):
        print(mu[ii], end="\n")


print()
Evec = np.array([(20, 30, 40), (50, 60, 70)], dtype=np.single)
mu = GetMu('water',Evec)
print(mu,type(mu))
