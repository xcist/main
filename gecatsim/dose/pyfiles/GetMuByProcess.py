# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

import numpy as np
from ctypes import *
import numpy.matlib
import gecatsim as xc
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.ReadMaterialFile import ReadMaterialFile
    
def GetMuByProcess(cfg,MaterialName = None,Evec = None,ProcessName = None,PairProductionFlag = None): 
    
    #print('Setting mu for for material %s...'%MaterialName)
    # It seems that for over 1024 kVp (1MV), the pair production flag is set.
    if PairProductionFlag is None:
        PairProductionFlag = np.max(Evec) > 1024
    
    # find the right directory TODO
    #if (MaterialDirectory is None) or len(MaterialDirectory)==0:
        #FIXME
        #SetupDirs
     #   pass
    MaterialDirectory = my_path.paths["material"] + "/"
    c_MaterialDirectory = c_char_p(bytes(MaterialDirectory, 'utf-8'))
    
    if MaterialName == 'unit':
        Mu = np.single(1.0)
    else:
        # Read the material file
        MaterialFile = my_path.find("material", MaterialName, '')
        Number, Density,AtomicNumbers,NormalizedMassFractions = ReadMaterialFile(MaterialFile)
        feval('C_Materials_CrossSectionDB_Initialize',cfg, c_MaterialDirectory, PairProductionFlag)
        # Get the normalized mass attenuation coefficients for each element in the material.
        #MAC = np.zeros((Evec.shape,Evec.shape))
        #MAC = np.matlib.repmat(MAC,4,1)
        MAC = feval('C_Materials_CrossSectionMAC_ByProc_Get',cfg, AtomicNumbers,NormalizedMassFractions,Evec)
        # Because the mass fractions passed are normalized by Materials_ElementalComposition_Get,
        # units are cm^2/g of material (not per gram of element).
        # Convert MAC to mu for each element in the material.
        Mu = MAC * Density
        # Units are cm^-1 (cm^2/g * g/cm^3)
        numE = int(np.asarray(Mu).size/ 4)# in matlab
        ProcessName = ProcessName.lower()
        if ProcessName == 'rayleigh' or ProcessName == 'r':
            ind = np.arange(numE)
        elif ProcessName == 'compton' or ProcessName =='c':
            ind = np.arange(numE,numE * 2)
        elif ProcessName == 'photoelectric' or ProcessName == 'pe':
            ind = np.arange(numE * 2,numE * 3)
        elif ProcessName == 'pairproduction' or ProcessName =='pp':
            ind = np.arange(numE * 3,numE * 4)
        else:
            ind = np.arange(numE * 2,numE * 3)
        Mu = Mu[ind]
        Mu = np.reshape(Mu, np.array(Evec).shape, order="F")
    
    #print('... done setting mu for for %s...'%MaterialName)
    return Mu
