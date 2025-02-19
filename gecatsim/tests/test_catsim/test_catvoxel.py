import unittest
import numpy as np
from unittest.mock import patch, MagicMock, mock_open
from io import StringIO
from gecatsim.pyfiles.catvoxel import catvoxel
from gecatsim.pyfiles import CommonTools

class TestCatVoxel(unittest.TestCase):

    @patch('gecatsim.pyfiles.catvoxel.my_path.find')
    @patch('gecatsim.pyfiles.catvoxel.json.load')
    @patch('gecatsim.pyfiles.catvoxel.GetMu')
    @patch('gecatsim.pyfiles.catvoxel.MakeAllVolumes')
    def test_catvoxel(self, mock_make_all_volumes, mock_get_mu, mock_json_load, mock_find):
        # Mocking the cfg object
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")

        cfg.phantom.filename = '../phantom/CatSimLogo_1024/CatSim_logo_1024.json'
        cfg.make_img_kv = 120
        cfg.recon_size = 256
        cfg.recon_planes = 128
        cfg.recon_fov = 500
        cfg.recon_slice_thickness = 1
        cfg.recon_xcenter = 0
        cfg.recon_ycenter = 0
        cfg.recon_zcenter = 0
        cfg.vol_os = 2
        cfg.material_volumes = True
        cfg.write_vp = True

        # Mocking the setMaterial method
        cfg.phantom.setMaterial = MagicMock()

        mock_find.return_value = 'phantom_file_path'
        mock_json_load.return_value = {
            'mat_name': ['Material1', 'Material2'],
            'n_materials': 2
        }
        mock_get_mu.side_effect = [np.array([0.1, 0.2]), np.array([0.3, 0.4])]
        mock_make_all_volumes.return_value = np.zeros((256, 256, 128, 2))

        # Mocking file operations
        with patch('builtins.open', mock_open(read_data='{"mat_name": ["Material1", "Material2"], "n_materials": 2}')):
            with patch('os.path.splitext', return_value=('phantom_file_path', '.json')):
                with patch('io.open', lambda *args, **kwargs: StringIO() if 'w' in args else open(*args, **kwargs)):
                    catvoxel('config_file', cfg)

        # Basic assertions to ensure the function runs without errors
        mock_find.assert_called_once()
        mock_json_load.assert_called_once()
        mock_get_mu.assert_any_call('Material1', 120)
        mock_get_mu.assert_any_call('Material2', 120)
        mock_make_all_volumes.assert_called_once()
        cfg.phantom.setMaterial.assert_called_once()

if __name__ == '__main__':
    unittest.main()