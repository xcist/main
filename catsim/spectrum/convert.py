import numpy as np
from scipy.io import loadmat

inppath = 'xcist_kVp80_tar7_bin1.mat'
inp = loadmat(inppath)
Eng = inp['newE'][0]
spec = inp['newI2'][0]

with open(inppath.replace('mat','dat'), 'w') as f:
    f.write('\n')
    f.write('%d\n'%len(Eng))
    for i in range(len(Eng)):
        f.write("%e,%e\n"%(Eng[i], spec[i]))
