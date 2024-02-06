from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.Phantom_Polygonal import (
    set_materials,
    C_Phantom_Polygonal_Clear
)
from unittest.mock import patch, MagicMock
from gecatsim.pyfiles.GetMu import GetMu
from ctypes import *
from numpy.ctypeslib import ndpointer

class Test_Phantom_Polygonal_set_materials(unittest.TestCase):
    def test_Phantom_Polygonal_set_materials(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")
        materialList = [
            'ncat_water',
            'ncat_muscle',
            'ncat_lung',
            'ncat_dry_spine',
            'ncat_dry_rib',
            'ncat_adipose',
            'ncat_blood',
            'ncat_heart',
            'ncat_kidney',
            'ncat_liver',
            'ncat_lymph',
            'ncat_pancreas',
            'ncat_intestine',
            'ncat_skull',
            'ncat_cartilage',
            'ncat_brain',
            'ncat_spleen',
            'ncat_iodine_blood',
            'ncat_iron',
            'ncat_pmma',
            'ncat_aluminum',
            'ncat_titanium',
            'ncat_air',
            'ncat_graphite',
            'ncat_lead',
            'ncat_breast_mammary',
            'ncat_skin',
            'ncat_iodine',
            'ncat_eye_lens',
            'ncat_ovary',
            'ncat_red_marrow',
            'ncat_yellow_marrow',
            'ncat_testis',
            'ncat_thyroid',
            'ncat_bladder']

        cfg.spec.Evec = range(10, 70, 10)
        # nMat = len(materialList)
        # Mus = np.zeros([len(cfg.spec.Evec), nMat], dtype=np.float64)
        # for i in range(nMat):
        #     Mus[:, i] = GetMu(materialList[i], cfg.spec.Evec) / 10

        cfg.clib.set_material_info = MagicMock()
        set_materials(cfg,materialList)

        assert cfg.clib.set_material_info.call_count == 1


class Test_C_Phantom_Polygonal_Clear(unittest.TestCase):
    def test_C_Phantom_Polygonal_Clear(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                    "../examples/cfg/Protocol_Sample_axial")

        num_polygons = 5
        cfg.clib.clear_polygonalized_phantom = MagicMock()
        C_Phantom_Polygonal_Clear(cfg, num_polygons)

        assert cfg.clib.clear_polygonalized_phantom.call_count == 1