from gecatsim.pyfiles.CommonTools import *
from gecatsim.reconstruction.pyfiles.createHSP import *
import unittest.mock

class Test_createHSP(unittest.TestCase):
    def test_createHSP_1(self):
        length = 4

        kerneltype = "R-L"
        fft_f = createHSP(length, kerneltype)

        assert(fft_f == [0.+0.04735763j, 0.+0.10671021j, 0.+0.25j,       0.+0.39328979j,
                            0.+0.45264237j, 0.+0.39328979j, 0.+0.25j,       0.+0.10671021j]).all

    def test_createHSP_2(self):
        length = 4

        kerneltype = "S-L"
        fft_f = createHSP(length, kerneltype)

        assert(fft_f == [0.+0.05403796j, 0.01350949+0.10711584j,  0.+0.21615186j,
                            -0.01350949+0.2981689j, 0.+0.32422779j, 0.01350949+0.2981689j,
                            0.+0.21615186j, -0.01350949+0.10711584j]).all

    def test_createHSP_3(self):
        length = 4

        kerneltype = "R-L"
        fft_f = createHSP(length, kerneltype)

        assert(fft_f == [0.+0.04735763j, 0.+0.10671021j, 0.+0.25j, 0.+0.39328979j,
                            0.+0.45264237j, 0.+0.39328979j, 0.+0.25j, 0.+0.10671021j]).all

    def test_createHSP_4(self):
        length = 4

        kerneltype = "soft"
        fft_f = createHSP(length, kerneltype)

        assert(fft_f == [0.+0.00037388j, 0.+0.10187408j, 0.+0.11392834j, 0.+0.06122712j,
                            0.+0.08161577j, 0.+0.17080719j, 0.+0.20344346j, 0.+0.12499887j]).all

    def test_createHSP_5(self):
        length = 4

        kerneltype = "standard"
        fft_f = createHSP(length, kerneltype)

        assert(fft_f == [0.+0.00037388j, 0.+0.11672395j, 0.+0.18574513j, 0.+0.16560513j,
                            0.+0.22075169j, 0.+0.27847859j, 0.+0.23309877j, 0.+0.12499887j]).all

    def test_createHSP_6(self):
        length = 4

        kerneltype = "bone"
        fft_f = createHSP(length, kerneltype)

        assert(fft_f == [0.+3.73875000e-04j, 0.+1.31061320e-01j, 0.+2.92059934e-01j,
                            0.+4.56658477e-01j, 0.+6.08725902e-01j, 0.+4.37871184e-01j,
                            0.+2.61730633e-01j, 0.+1.24998875e-01j]).all

    def test_createHSP_7(self):
        length = 4

        kerneltype = "fakekernel"
        try:
            fft_f = createHSP(length, kerneltype)
        except Exception as ex:
            exception_string = str(ex)
            assert(exception_string=='******** Error! An unsupported kernel was specified: fakekernel. ********')

