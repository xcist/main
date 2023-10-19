from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.randpf import randpf
import unittest.mock

class Test_randpf(unittest.TestCase):
    def test_randpf1(self):
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

        check_value(out_array1)
    def test_randpf2(self):
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
        check_value(out_array2)

    def test_randpf3(self):
        input_array3 = np.float32([[77.60277,   70.674126,  28.36824,   59.245857 ],
                        [18.446276,  18.401636,  45.93234,   20.344593 ],
                        [84.89199,   69.48412,   24.458265,  55.101936 ],
                        [ 3.5849595, 24.861238,   5.987396,  13.479603 ],
                        [86.30522,   4.846194,  74.6355,    33.699875 ]]
                                  )
       # input_array3 = np.float32(input_array3)
        input_array3[1:4, 1:3] = 0
        out_array3 = randpf(input_array3)
        assert (out_array3 == [[73., 80., 22., 55.],
                               [17.,  0.,  0., 24.],
                               [89.,  0.,  0., 64.],
                               [ 4.,  0.,  0., 11.],
                               [91.,  4., 76., 48.]]

                ).any()
        check_value(out_array3)
