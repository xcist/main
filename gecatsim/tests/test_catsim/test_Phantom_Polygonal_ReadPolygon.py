from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.Phantom_Polygonal_ReadPolygon import extract_polygonal_objects

class Test_Phantom_Polygonal_ReadPolygon(unittest.TestCase):
    def test_Phantom_Polygonal_ReadPolygon(self):
        file_path = "../../phantom/reduced_female_10yr_lung_lesions.nrb"
        _objects = extract_polygonal_objects(file_path)
        assert _objects is not None
        assert _objects['vertices'] is not None
        assert  _objects['num_triangles'] is not None
        assert _objects['materialId'] is not None



if __name__ == '__main__':
    unittest.main()
