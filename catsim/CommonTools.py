# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import numpy as np
import numpy.matlib as nm
import ctypes, struct
import os

'''
Common tool functions.
Mingye Wu, GE Research

'''

def check_value(a):
    #return
    print(type(a).__name__)
    if type(a) is np.ndarray:
        print(a.shape, a.dtype)
    print(a,'\n')
    
def make_col(a):
    a = a.reshape(a.size, 1)
    return a
        
def feval(funcName, *args):
    try:
        md = __import__(funcName)
    except:
        md = __import__("catsim."+funcName, fromlist=[funcName])  # equal to: from catsim.foo import foo
    return eval("md."+funcName)(*args)

def get_path():
    # Locate paths of lib and data.
    myPath = emptyCFG()
    myPath.main = os.path.dirname(os.path.abspath(__file__))
    myPath.cfg = myPath.main+'/cfg'
    myPath.lib = myPath.main+'/lib'
    myPath.bowtie = myPath.main+'/data/bowtie'
    myPath.material = myPath.main+'/data/material'
    myPath.phantom = myPath.main+'/data/phantom'
    myPath.spectrum = myPath.main+'/data/spectrum'
    return myPath

def load_C_lib():
    myPath = get_path()

    # add lib path to environment value "PATH" for depending DLLs
    os.environ["PATH"] = myPath.lib+';'+os.environ["PATH"]

    # load C/C++ lib
    ll = ctypes.cdll.LoadLibrary
    if os.name == "nt":
        libFile = "libcatsim64.dll"
    else:
        libFile = "libcatsim.so"
    clib = ll(myPath.lib + "/" + libFile)
    
    return clib

class emptyCFG:
    pass
    
class CFG:
    def __init__(self, *para):
        # initialize cfg: defaults, paths, and C lib
        cfg = source_cfg("Phantom_Default")
        cfg = source_cfg("Scanner_Default", cfg)
        cfg = source_cfg("Protocol_Default", cfg)
        cfg = source_cfg("Physics_Default", cfg)
        cfg = source_cfg("Recon_Default", cfg)
        cfg.resultsName = "simulation_test"

        cfg.path = get_path()
        cfg.clib = load_C_lib()
        
        # source cfgFiles if para are defined
        # note: the later cfgFile overrides the former ones
        for cfgFile in para:
            cfg = source_cfg(cfgFile, cfg)
        
        self.pass_cfg_to_self(cfg)

    def pass_cfg_to_self(self, cfg):
        # add or override cfg attributes to self
        for name1, value1 in vars(cfg).items():
            if not hasattr(self, name1):
                setattr(self, name1, value1)
            else:
                for name2, value2 in eval("vars(cfg.%s).items()" % name1):
                    setattr(getattr(self, name1), name2, value2)
            
    def load(self, cfgFile):
        cfg = source_cfg(cfgFile)
        self.pass_cfg_to_self(cfg)
    
def source_cfg(*para):
    '''
    First para must be cfg filename.
    Second para is optional, if defined and is cfg, attr will be added to cfg.
    Calling source_cfg(cfgFile, cfg) will add or override attributes to cfg.
    '''
    # find cfg file
    cfgFile = para[0]
    if not os.path.isfile(cfgFile):
        cfgPath = get_path().cfg + "/"
        if os.path.isfile(cfgFile + ".cfg"):
            cfgFile += ".cfg"
        elif os.path.isfile(cfgPath + cfgFile):
            cfgFile = cfgPath + cfgFile
        elif os.path.isfile(cfgPath + cfgFile + ".cfg"):
            cfgFile = cfgPath + cfgFile + ".cfg"
        else:
            raise Exception("Cannot find %s or %s.cfg" %(cfgFile, cfgFile))
    
    # cfg is initialized before sourcing cfgFile
    if len(para)<2:
        cfg = emptyCFG()
    else:
        cfg = para[1]
        
    # initialize structs in cfg and structs
    attrList = ['sim', 'det', 'detNew', 'src', 'srcNew', 'spec', 'protocol', 'scanner', 'phantom', 'physics', 'recon']
    for attr in attrList:
        if not hasattr(cfg, attr):
            setattr(cfg, attr, emptyCFG())
        if not attr in dir():
            exec("%s = emptyCFG()" % attr)
        
    # execute scripts in cfg file
    exec(open(cfgFile).read())
    
    # add or override the attributes in the original cfg
    for attr in attrList:
        for name, value in eval("vars(%s).items()" % attr):
            setattr(getattr(cfg, attr), name, value)
    
    return cfg
    
def vectornorm(xyz):
    if xyz.shape[0]!=3:
        print('ERROR : argument of vectornorm has to be of size 3 x n');
        return
    else:
        norms = np.sqrt(np.square(xyz).sum(axis=0))
        norms = make_col(norms)
        return norms
        
def overlap(x0, y0, x1):
    # length
    n0 = len(x0)
    n1 = len(x1)
    
    # boundaries
    b0 = get_vector_boundaries(x0)
    b1 = get_vector_boundaries(x1)
    
    # default
    y1 = np.zeros(x1.shape)
    
    # pre-loop, find the start boundaries
    i = 0
    j = 0
    previous = b1[j]
    while b0[i] < previous:
        i += 1
        if i>=n0:
            return y1
    if i == 0:
        previous = b0[0]
        while b1[j] < previous:
            j += 1
            if j >= n1:
                return y1    
    
    # main loop
    while j < n1:
        if b0[i+1] < b1[j+1]:
            y1[j] += y0[i]*(b0[i+1]-previous)/(b1[j+1]-b1[j])
            previous = b0[i+1]
            i += 1
            if i >= n0:
                return y1
        else:
            y1[j] += y0[i]*(b1[j+1]-previous)/(b1[j+1]-b1[j])
            previous = b1[j+1]
            j += 1  
    
    return y1

def get_vector_boundaries(x):
    # x can be scalar, vector, or [n, 1] array
    
    if len(x) == 1:
        b = np.array([x*(1-1e-6), x*(1+1e-6)])
    else:
        b = (x[0:-1]+x[1:])/2
        b = np.concatenate(([x[0]-0.5*(x[1]-x[0])], b, [x[-1]+0.5*(x[-1]-x[-2])]))
    return b

def rawread(fname, dataShape, dataType):
    # dataType is for numpy, can be 'float' or 'single', 'double', 'int', 'uint' ONLY
    #          they are single, double, int32, uin32
    with open(fname, 'rb') as fin:
        data = fin.read()
    
    switcher = {'float': ['f', 4, np.single], 'single': ['f', 4, np.single], 'double': ['d', 8, np.double], 'int': ['i', 4, np.int32], 'uint': ['I', 4, np.uint32]}
    fmt = switcher[dataType]
    data = struct.unpack("%d%s" % (len(data)/fmt[1], fmt[0]), data)
    
    data = np.array(data, dtype=fmt[2])
    if dataShape:
        data = data.reshape(dataShape)
    
    return data

def rawwrite(fname, data):
    with open(fname, 'wb') as fout:
        fout.write(data)