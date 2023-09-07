from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.randpf import randpf
import unittest.mock

class test_randpf(unittest.TestCase):
    def Test_randpf(self):
        input_array1 = np.float32([[30., 28., 10., 4.],
                           [11., 0., 0., 37.],
                           [79., 0., 0., 82.],
                           [47., 0., 0., 74.],
                           [73., 8., 69., 91.]]
                        )

        input_array1[1:4, 1:3] = 0
        out_array1 = randpf(input_array1)
        assert (out_array1 == np.array(
            [[29., 17., 15., 2.],
             [15., 0., 0., 38.],
             [64., 0., 0., 58.],
             [61., 0., 0., 74.],
             [63., 6., 72., 100.]]
                                )).any()

    def test2_randpf(self):
        input_array2 = np.float32([[12.244902, 43.070995, 1.0779757, 63.588734],
                                       [20.174131, 0., 0., 93.01599],
                                       [75.9809, 0., 0., 28.409391],
                                       [57.16415, 0., 0., 21.542444],
                                       [5.9953947, 37.528004, 2.6598074, 57.94161]]
                                      )



        input_array2[1:4, 1:3] = 0
        out_array2 = randpf(input_array2)
        assert (out_array2 == np.array(
            [[11., 43.,  0., 59.],
             [24.,  0.,  0., 74.],
             [80.,  0.,  0., 28.],
             [47.,  0.,  0., 17.],
             [5., 35.,  5., 55.]]
        )).any()

