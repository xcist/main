import unittest.mock
from unittest.mock import patch

import numpy as np

from gecatsim.pyfiles.CommonTools import load_C_lib
from gecatsim.pyfiles.GetMu import GetMu


class TestGetMu(unittest.TestCase):

    def test_get_mu_precalculated_values(self):
        Evec = range(10, 70, 10

        mu = GetMu('water', Evec)

        assert mu == [5.325486183166504, 0.8101319670677185, 0.3757021129131317, 0.26842963695526123,
                      0.22698049247264862, 0.20569877326488495]

        Evec = np.array(
            [3.25, 9.25, 15.25, 21.25, 27.25, 33.25, 39.25, 45.25, 51.25, 57.25, 63.25, 69.25, 75.25, 81.25, 87.25,
             93.25,
             99.25, 105.25, 111.25, 117.25])

        mu = GetMu('Al', Evec)

        assert (mu == np.array(
            [1.7148381e+03, 8.8936996e+01, 2.0453079e+01, 7.8083868e+00, 3.9207914e+00, 2.3513255e+00, 1.6005120e+00,
             1.1929718e+00, 9.5403385e-01, 8.0138618e-01, 7.0024604e-01, 6.2837154e-01, 5.7653737e-01, 5.3749001e-01,
             5.0646776e-01, 4.8200923e-01, 4.6231583e-01, 4.4530702e-01, 4.3088070e-01, 4.1860345e-01],
            dtype=np.single)).all()

    @patch('gecatsim.pyfiles.GetMu.load_C_lib', create=True)
    def test_get_mu_calls_load_clib(self, load_C_lib_mock):
        clib = load_C_lib()
        load_C_lib_mock.return_value = clib
        Evec = range(10, 70, 10)

        mu = GetMu('water', Evec)

        assert load_C_lib_mock.call_count == 1
