import unittest
import numpy as np
from unittest.mock import patch, mock_open
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.catsimfilter import catsimfilter

class TestCatsimfilter(unittest.TestCase):

    @patch('os.path.exists')
    @patch('os.path.isdir')
    @patch('os.path.isfile')
    @patch('os.listdir')
    @patch('gecatsim.pyfiles.Spectrum.spectrum_read')
    @patch('gecatsim.pyfiles.GetMu.ReadMaterialFile')
    @patch('gecatsim.pyfiles.GetMu.GetMu')
    @patch('builtins.open', new_callable=mock_open, read_data="3\n10,1\n20,2\n30,3\n")
    def test_catsimfilter(self, mock_open_fn, mock_get_mu, mock_read_material_file,
                          mock_spectrum_read, mock_listdir, mock_isfile, mock_isdir, mock_exists):

        # Setup cfg
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

        # Mock filesystem and function behavior
        mock_exists.side_effect = lambda path: path in [
            'tungsten_tar7_120_unfilt.dat',
            'gecatsim/spectrum/tungsten_tar7_120_unfilt.dat'
        ]
        mock_isdir.return_value = False
        mock_isfile.return_value = True
        mock_spectrum_read.return_value = (np.array([10, 20, 30]), np.array([1, 2, 3]), 0)
        mock_read_material_file.return_value = (2, 1.0, [1, 2], [0.5, 0.5])
        mock_get_mu.return_value = np.array([0.1, 0.2, 0.3])

        catsimfilter(cfg)

        assert mock_open_fn.called
        handle = mock_open_fn()
        written = [call.args[0] for call in handle.write.call_args_list]

        assert '3\n' in written
        assert any('10.0' in line for line in written)
        assert any('20.0' in line for line in written)
        assert any('30.0' in line for line in written)

if __name__ == '__main__':
    unittest.main()