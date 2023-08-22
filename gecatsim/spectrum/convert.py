import numpy as np
import re
from scipy.io import loadmat
from glob import glob
import matplotlib.pyplot as plt

allpath =  glob("*mat")
for inppath in allpath:
    inp = loadmat(inppath)
    Eng = inp['newE'][0]
    spec = inp['newI'][0]
    plt.plot(Eng, spec);plt.show()
    
    with open(inppath.replace('mat','dat'), 'w') as f:
        f.write('\n')
        f.write('%d\n'%len(Eng))
        for i in range(len(Eng)):
            f.write("%e,%e\n"%(Eng[i], spec[i]))
        # write angle
        ang = float(re.findall("\d+\.*\d*", inppath.split('_')[2])[0])
        f.write("%e"%ang)
