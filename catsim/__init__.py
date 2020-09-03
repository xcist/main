__version__ = '0.1.0'
from .CatSim import CatSim
from .GetMu import GetMu
from .CommonTools import rawread, rawwrite, check_value, get_path, CFG, source_cfg

def help():
    print('''e.g. with "import catsim as xc", the following functions/classes can be called:
    xc.CatSim(cfgFile0, cfgFile1, cfgFile2, ...)
    xc.GetMu(materialFile, Evec)
    xc.rawread(fname, dataShape, dataType)
    xc.rawwrite(fname, data)
    xc.check_value(var)
    xc.get_path()
    xc.CFG(cfgFile0, cfgFile1, cfgFile2, ...)
    xc.source_cfg(cfgFile[, cfg])
    ''')