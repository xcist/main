from gecatsim.pyfiles.C_Materials_CrossSectionMAC_ByProc_Get import C_Materials_CrossSectionMAC_ByProc_Get
import unittest.mock
import gecatsim.pyfiles.CommonTools as c
from unittest.mock import patch, MagicMock

class TestC_Materials_CrossSectionMAC_ByProc_Get(unittest.TestCase):

    def test_C_Materials_CrossSectionMAC_ByProc_Get_calls_GetCrossSectionByProcessMAC(self):
        cfg = c.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic", "../examples/cfg/Protocol_Sample_axial")
        cfg.clib.GetCrossSectionByProcessMAC = MagicMock()
        C_Materials_CrossSectionMAC_ByProc_Get(cfg,[1,2],[1.0,2.0],[1.0,2.0],[1.0,2.0])
        assert cfg.clib.GetCrossSectionByProcessMAC.call_count == 1
