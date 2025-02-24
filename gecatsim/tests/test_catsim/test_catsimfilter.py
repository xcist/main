import unittest
import numpy as np
from unittest.mock import patch, mock_open, MagicMock, call
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.Spectrum import spectrum_read
from gecatsim.pyfiles.catsimfilter import catsimfilter

class TestCatsimfilter(unittest.TestCase):

    def test_catsimfilter(self):
        cfg = CFG("../examples/cfg/Phantom_Sample",
                  "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")
        cfg.flat_filters = ['filter1', 'filter2']
        cfg.spectrum_filename = 'tungsten_tar7_120_unfilt.dat'
        cfg.spectrum_filtered_string = '_test'
        cfg.det.gammas = np.array([0.1, 0.2, 0.3])
        cfg.det.alphas = np.array([0.1, 0.2, 0.3])
        cfg.det.totalNumCells = 3
        cfg.spec.nEbin = 3
        cfg.protocol.flatFilter = ['material1', 0.1, 'material2', 0.2]
        cfg.src.filterTrans = np.ones((3, 3))

        @patch('os.path.exists')
        @patch('os.path.isdir')
        @patch('os.path.isfile')
        @patch('os.listdir')
        @patch('gecatsim.pyfiles.Spectrum.spectrum_read')
        @patch('gecatsim.pyfiles.GetMu.ReadMaterialFile')
        @patch('gecatsim.pyfiles.GetMu.GetMu')
        @patch('builtins.open', new_callable=mock_open, read_data="3\n10,1\n20,2\n30,3\n")
        def run_test(mock_open, mock_get_mu, mock_read_material_file, mock_spectrum_read, mock_listdir, mock_isfile, mock_isdir, mock_exists):
            # Mock the os.path.exists to return True for the spectrum file
            def exists_side_effect(path):
                if path in ['tungsten_tar7_120_unfilt.dat', 'gecatsim/spectrum/tungsten_tar7_120_unfilt.dat']:
                    return True
                return False
            mock_exists.side_effect = exists_side_effect

            mock_isdir.return_value = False
            mock_isfile.return_value = True

            mock_spectrum_read.return_value = (np.array([10, 20, 30]), np.array([1, 2, 3]), 0)
            mock_read_material_file.return_value = (2, 1.0, [1, 2], [0.5, 0.5])
            mock_get_mu.return_value = np.array([0.1, 0.2, 0.3])

            catsimfilter(cfg)

            mock_open.assert_any_call('filtered\\tungsten_tar7_120_unfilt.dat', 'w')

            # Check if the file was written correctly
            expected_calls = [
                call('tungsten_tar7_120_unfilt.dat', 'r'),
                call().__iter__(),
                call('filtered\\tungsten_tar7_120_unfilt.dat', 'w'),
                call().__enter__(),
                call().write('3\n'),
                call().write('10.0,[0.37138876 0.7420587  1.1126724 ]\n'),
                call().write('20.0,[0.3714966  0.74225223 1.1129498 ]\n'),
                call().write('30.0,[0.37168625 0.7425928  1.1134382 ]\n'),
                call().__exit__(None, None, None)
            ]
            mock_open.assert_has_calls(expected_calls, any_order=True)

        run_test()

if __name__ == '__main__':
    unittest.main()