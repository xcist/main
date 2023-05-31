import numpy as np
from ctypes import *
from numpy.ctypeslib import ndpointer
    
def C_Materials_CrossSectionMAC_ByProc_Get(cfg, AtomicNumbers = None,MassFractions = None,Energies = None,MAC = None): 
    NumberOfElements = np.array(AtomicNumbers).size
    #MassFractions = np.single(MassFractions)
    #AtomicNumbers = np.int32(AtomicNumbers)
    AtomicNumbers = (c_int*NumberOfElements)(*AtomicNumbers)
    MassFractions = (c_float*NumberOfElements)(*MassFractions)
    neng = len(Energies)
    Energies = (c_float*neng)(*Energies)
    #Energies = np.single(Energies)
    #MassFractions = np.single(Mac)
    NumberOfEnergies = np.asarray(Energies).size
    #print('Getting the MACs for %d elements from C code.'%NumberOfElements)
    # Get the mass attenuation coefficients for each element in the list of all four processes.
    #cfgnew = cfg.get_current_cfg()
    fun = cfg.clib.GetCrossSectionByProcessMAC
    fun.argtypes = [c_int, POINTER(c_int), POINTER(c_float), c_int, POINTER(c_float), POINTER(c_float)]
    MAC = (c_float*(4*NumberOfEnergies))()
    fun.restype = None
    # Load the cross-section database into the simulator (as members of the 'CrossSection' class).
    fun(NumberOfElements,AtomicNumbers,MassFractions,NumberOfEnergies,Energies,MAC)
    return np.array(MAC)
