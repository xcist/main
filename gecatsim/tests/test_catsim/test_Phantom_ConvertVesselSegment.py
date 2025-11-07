import unittest
from gecatsim.pyfiles.Phantom_ConvertVesselSegment import objs, endface

class TestVesselSegmentModeling(unittest.TestCase):

    def test_clip_structure(cls):
        """Ensure each clip matrix has shape (4, 2)."""
        for clip in objs['clip']:
            cls.assertEqual(clip.shape, (4, 2))

    def test_endface_keys(cls):
        """Check that each endface entry contains required keys."""
        expected_keys = {'center', 'R', 'curv_vec', 'cent_line_vec'}
        for entry in endface.values():
            cls.assertTrue(expected_keys.issubset(entry.keys()))

    def test_params_structure(cls):
        """Verify that each params entry has 15 elements."""
        for params in objs['params']:
            cls.assertEqual(len(params), 15)

if __name__ == '__main__':
    unittest.main()
