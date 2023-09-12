from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.ReadMaterialFile import ReadMaterialFile
import unittest.mock

class test_ReadMaterialFile(unittest.TestCase):
    def Test_ReadMaterialFile(self):
        (numberOfElements, density, atomicNumbers, massFractions) = ReadMaterialFile('../material/water')
        print(numberOfElements, density, atomicNumbers, massFractions)

        assert numberOfElements == 2
        assert density == 1.0
        assert (atomicNumbers == [1, 8]).all()
        assert massFractions == [0.111902, 0.888098].all()
