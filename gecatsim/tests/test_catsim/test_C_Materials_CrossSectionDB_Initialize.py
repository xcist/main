from gecatsim.pyfiles.C_Materials_CrossSectionDB_Initialize import C_Materials_CrossSectionDB_Initialize
import unittest.mock
import gecatsim.pyfiles.CommonTools as c
from unittest.mock import patch, MagicMock

class TestC_Materials_CrossSectionDB_Initialize(unittest.TestCase):

    def test_C_Materials_CrossSectionDB_Initialize_calls_clib_InitializeCrossSectionDB(self):
        cfg = c.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic", "../examples/cfg/Protocol_Sample_axial")
        cfg.clib.InitializeCrossSectionDB = MagicMock()
        C_Materials_CrossSectionDB_Initialize(cfg)
        assert cfg.clib.InitializeCrossSectionDB.call_count == 1
