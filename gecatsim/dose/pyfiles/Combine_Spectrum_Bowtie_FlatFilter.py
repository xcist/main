import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
    
def Combine_Spectrum_Bowtie_FlatFilter(cfg): 
    print('Combining spectrum with bowtie transmission and flat filter transmission...')
    # The spectrum is in units of photons/view/mm^2 per energy bin at 1m.
    if cfg.spec.Ivec.shape == cfg.src.filterTrans.shape:
        cfg.spec.netIvec = np.multiply(cfg.spec.Ivec, cfg.src.filterTrans)
    else:
        raise Exception('spec.Ivec, bowtie.transVec and FiltrationTransVec(flat filter) should have the same dimension')
    
    print('... done combining spectrum with bowtie transmission and flat filter transmission.')
    
    return cfg
