import unittest.mock
from unittest.mock import MagicMock

import gecatsim.pyfiles.CommonTools as c
from gecatsim.pyfiles.C_Materials_CrossSectionDB_Initialize import C_Materials_CrossSectionDB_Initialize


class TestC_Materials_CrossSectionDB_Initialize(unittest.TestCase):

    def test_c_materials_cross_sectiondb_initialize_calls_clib_initializecross_sectiondb(self):
        cfg = c.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                    "../examples/cfg/Protocol_Sample_axial")

        cfg.clib.InitializeCrossSectionDB = MagicMock()
        C_Materials_CrossSectionDB_Initialize(cfg)

        assert cfg.clib.InitializeCrossSectionDB.call_count == 1
