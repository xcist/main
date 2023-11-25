from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.ReadMaterialFile import ReadMaterialFile
from gecatsim.material import *
import unittest.mock

class Test_ReadMaterialFile(unittest.TestCase):
    def test_ReadMaterialFile(self):

        # print("/n")
        # print(os.path.exists('../../material/water'))

        (numberOfElements, density, atomicNumbers, massFractions) = ReadMaterialFile('gecatsim/material/water')

        assert numberOfElements == 2
        assert density == 1.0
        assert atomicNumbers == [1, 8]
        assert massFractions == [0.111902, 0.888098]
