from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.ReadMaterialFile import ReadMaterialFile
import unittest.mock

class Test_ReadMaterialFile(unittest.TestCase):
    def test_ReadMaterialFile(self):
        (numberOfElements, density, atomicNumbers, massFractions) = ReadMaterialFile('water')

        assert numberOfElements == 2
        assert density == 1.0
        assert atomicNumbers == [1, 8]
        assert massFractions == [0.111902, 0.888098]
