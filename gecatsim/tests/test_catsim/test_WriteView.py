import unittest
import os
import numpy as np
from gecatsim.pyfiles.WriteView import WriteView

class TestWriteView(unittest.TestCase):

    def setUp(self):
        self.filename = 'test_view.dat'
        self.view_data = np.array([1.0, 2.0, 3.0], dtype='float32')

    def tearDown(self):
        if os.path.exists(self.filename):
            os.remove(self.filename)
        for i in range(1, 4):
            if os.path.exists(f"{self.filename}.{i}"):
                os.remove(f"{self.filename}.{i}")

    def test_write_view_single_file(self):
        WriteView(self.filename, 1, self.view_data, SplitSinogramPerView=0)
        self.assertTrue(os.path.exists(self.filename))
        with open(self.filename, 'rb') as f:
            data = np.fromfile(f, dtype='float32')
            np.testing.assert_array_equal(data, self.view_data)

    def test_write_view_split_files(self):
        for i in range(1, 4):
            WriteView(self.filename, i, self.view_data, SplitSinogramPerView=1)
            self.assertTrue(os.path.exists(f"{self.filename}.{i}"))
            with open(f"{self.filename}.{i}", 'rb') as f:
                data = np.fromfile(f, dtype='float32')
                np.testing.assert_array_equal(data, self.view_data)

    def test_write_view_append(self):
        WriteView(self.filename, 1, self.view_data, SplitSinogramPerView=0)
        WriteView(self.filename, 2, self.view_data, SplitSinogramPerView=0)
        self.assertTrue(os.path.exists(self.filename))
        with open(self.filename, 'rb') as f:
            data = np.fromfile(f, dtype='float32')
            np.testing.assert_array_equal(data, np.concatenate((self.view_data, self.view_data)))

if __name__ == '__main__':
    unittest.main()