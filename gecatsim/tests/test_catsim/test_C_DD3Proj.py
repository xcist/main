import unittest.mock
from unittest.mock import patch

import numpy as np

from gecatsim.pyfiles.C_DD3Proj import DD3Proj
from gecatsim.pyfiles.CommonTools import load_C_lib


class TestC_DD3Proj(unittest.TestCase):
    # --------- Set parameters
    # source coordinates
    x0 = 0
    y0 = 550
    z0 = 0

    # detector dimension and center coordinates
    nrdetcols = 200
    nrdetrows = 16

    xds = np.linspace(-nrdetcols / 2 + 0.5, nrdetcols / 2 - 0.5, nrdetcols, dtype=np.single)
    yds = -400 * np.ones(nrdetcols, dtype=np.single)
    zds = np.linspace(-nrdetrows / 2 + 0.5, nrdetrows / 2 - 0.5, nrdetrows, dtype=np.single)

    # original image and view settings
    dzdx = 1

    imgXoffset = 0
    imgYoffset = 0
    imgZoffset = 0

    nrviews = 315
    viewangles = np.linspace(0, np.pi * 2, nrviews, dtype=np.single)  # view angle of each view
    zshifts = np.zeros(nrviews, dtype=np.single)  # z-position of each view

    nrcols = 51
    nrrows = 51
    nrplanes = 30

    pOrig = np.zeros((nrrows, nrcols, nrplanes), dtype=np.single)
    pOrig[25, 25, :] = 1
    pOrig[25, 0, :] = 1
    pOrig[0, 0, :] = 1

    def test_c_dd3proj_precalculated_values(self):
        # --------- Run DD3Proj
        sinogram = DD3Proj(self.x0, self.y0, self.z0,
                           self.nrdetcols, self.nrdetrows,
                           self.xds, self.yds, self.zds,
                           self.dzdx,
                           self.imgXoffset, self.imgYoffset, self.imgZoffset,
                           self.viewangles,
                           self.zshifts,
                           self.nrviews,
                           self.nrcols, self.nrrows, self.nrplanes,
                           self.pOrig)

        # check values only in first dimension
        assert (sinogram[0][53] == [0.1430355, 0.14303426, 0.14303343, 0.1430326, 0.14303194, 0.14303152, 0.1430312,
                                    0.14303105, 0.14303103, 0.1430312, 0.14303152, 0.14303194, 0.1430326, 0.14303343,
                                    0.14303426, 0.1430355]).all
        assert (sinogram[0][54] == [1.0011778, 1.0011692, 1.0011634, 1.0011573, 1.001153, 1.00115, 1.0011476, 1.0011466,
                                    1.0011464, 1.0011476, 1.00115, 1.0011529, 1.0011573, 1.0011634, 1.0011692,
                                    1.0011778]).all
        assert (sinogram[0][55] == [0.7129221, 0.712916, 0.7129119, 0.7129076, 0.7129045, 0.71290225, 0.71290064,
                                    0.71289974, 0.7128998, 0.7129007, 0.71290225, 0.7129045, 0.7129076, 0.71291184,
                                    0.712916, 0.71292216]).all
        assert (sinogram[0][56] == [1.0010792, 1.0010707, 1.0010649, 1.0010589, 1.0010545, 1.0010514, 1.001049,
                                    1.001048, 1.001048, 1.001049, 1.0010514, 1.0010545, 1.0010589, 1.0010649, 1.0010707,
                                    1.0010792]).all
        assert (sinogram[0][57] == [0.6825226, 0.6825168, 0.68251276, 0.68250877, 0.6825057, 0.68250364, 0.68250203,
                                    0.6825013, 0.6825013, 0.68250203, 0.68250364, 0.6825057, 0.68250877, 0.68251276,
                                    0.6825168, 0.6825226]).all

        assert (sinogram[0][99] == [0.86366373, 0.8636562, 0.86365116, 0.8636461, 0.8636423, 0.86363965, 0.8636376,
                                    0.86363673, 0.86363673, 0.8636376, 0.86363965, 0.8636423, 0.8636461, 0.86365116,
                                    0.8636562, 0.86366373]).all
        assert (sinogram[0][100] == [0.86366373, 0.8636562, 0.86365116, 0.8636461, 0.8636423, 0.86363965, 0.8636376,
                                     0.86363673, 0.86363673, 0.8636376, 0.86363965, 0.8636423, 0.8636461, 0.86365116,
                                     0.8636562, 0.86366373]).all

        for i in range(0, self.nrviews):
            assert (sinogram[i][99] != np.zeros([self.nrdetrows], dtype=np.single)).all
            assert (sinogram[i][100] != np.zeros([self.nrdetrows], dtype=np.single)).all

    @patch('gecatsim.pyfiles.C_DD3Proj.load_C_lib', create=True)
    def test_c_dd3proj_calls_load_clib(self, load_C_lib_mock):
        clib = load_C_lib()
        load_C_lib_mock.return_value = clib

        sinogram = DD3Proj(self.x0, self.y0, self.z0,
                           self.nrdetcols, self.nrdetrows,
                           self.xds, self.yds, self.zds,
                           self.dzdx,
                           self.imgXoffset, self.imgYoffset, self.imgZoffset,
                           self.viewangles,
                           self.zshifts,
                           self.nrviews,
                           self.nrcols, self.nrrows, self.nrplanes,
                           self.pOrig)

        assert load_C_lib_mock.call_count == 1
