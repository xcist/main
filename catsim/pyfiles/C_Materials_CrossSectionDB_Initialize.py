import numpy as np
from ctypes import *
from numpy.ctypeslib import ndpointer
    
def C_Materials_CrossSectionDB_Initialize(cfg, MaterialDirectory = None,PairProductionFlag = None): 
    #print('Loading material cross-sections from %s',MaterialDirectory ,'and setting the cross-section database in C global variables.')
    if PairProductionFlag:
        print('WARNING! Pair production flag is on.')
    
    #cfgnew = cfg.get_current_cfg()
    fun = cfg.clib.InitializeCrossSectionDB
    fun.argtypes = [POINTER(c_char), c_int]
    fun.restype = None
    # Load the cross-section database into the simulator (as members of the 'CrossSection' class).
    fun(MaterialDirectory,PairProductionFlag)

    return
    
