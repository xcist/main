# TODO: should this version match setup.py?
__version__ = '0.1.8'

from gecatsim.pyfiles.CatSim import CatSim
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import rawread, rawwrite, check_value, CFG, source_cfg


def help():
    print('''e.g. with "import catsim as xc", the following functions/classes can be called:
    xc.help()
    xc.CommonTools.my_path.extra_search_paths.append("a-path-with-data-files")
    xc.CatSim(cfgFile0, cfgFile1, cfgFile2, ...)
    xc.GetMu(materialFile, Evec)
    xc.rawread(fname, dataShape, dataType)
    xc.rawwrite(fname, data)
    xc.check_value(var)
    xc.CFG(cfgFile0, cfgFile1, cfgFile2, ...)
    xc.source_cfg(cfgFile[, cfg])
    ''')
