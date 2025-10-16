import unittest
from unittest.mock import patch, MagicMock, mock_open
import numpy as np
import tempfile
import os
import json
from gecatsim.pyfiles.catvoxel import catvoxel
from gecatsim.pyfiles.catvoxel import MakeAllVolumes
from gecatsim.pyfiles.catvoxel import PrintPhantomSetup
from gecatsim.pyfiles.catvoxel import dumpjson

class TestCatVoxel(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.test_file = os.path.join(self.test_dir.name, "phantom.ph")

    def tearDown(self):
        self.test_dir.cleanup()

    
    @patch('gecatsim.pyfiles.catvoxel.MakeAllVolumes')
    @patch('gecatsim.pyfiles.catvoxel.feval')
    def test_materials_and_volume_shape(self, mock_feval, mock_makeallvolumes):

        # Setup mock cfg
        cfg = MagicMock()
        cfg.material_volumes = True
        cfg.write_vp = False
        cfg.make_img_kv = [100]
        cfg.spec = MagicMock()
        cfg.phantom.callback = 'phantom_callback'
        cfg.phantom.Materials = ['Water', 'Bone', 'Air']
        cfg.phantom.numberOfMaterials = 3
        cfg.phantom.projectorCallback = 'C_Projector_Analytic'
        cfg.phantom.filename = 'phantom.ph'
        cfg.recon.imageSize = 64
        cfg.recon.sliceCount = 32
        cfg.recon.fov = 256.0
        cfg.recon.sliceThickness = 1.0
        cfg.recon.centerOffset = [0.0, 0.0, 0.0]
        cfg.vol_os = 1

        # Mock feval to return cfg unchanged
        mock_feval.side_effect = lambda cb, c: c

        # Expected shape: (3, 32, 64, 64)
        dummy_volume = np.ones((3, 32, 64,64), dtype=np.single)
        mock_makeallvolumes.return_value = dummy_volume

        # Run the function
        catvoxel(cfg)

        # Assertions
        self.assertEqual(cfg.phantom.Materials, ['Water', 'Bone', 'Air'])
        self.assertEqual(cfg.phantom.numberOfMaterials, 3)
        mock_makeallvolumes.assert_called_once()
        args, kwargs = mock_makeallvolumes.call_args
        volume_arg = args[0]
        self.assertEqual(volume_arg.shape, (3, 32, 64, 64))


    @patch('gecatsim.pyfiles.catvoxel.feval', return_value=MagicMock())
    @patch('gecatsim.pyfiles.catvoxel.GetMu', return_value=MagicMock())
    @patch('gecatsim.pyfiles.catvoxel.MakeAllVolumes')
    @patch('gecatsim.pyfiles.catvoxel.PrintPhantomSetup')
    @patch('gecatsim.pyfiles.catvoxel.dumpjson')
    @patch('gecatsim.pyfiles.catvoxel.xc.rawwrite')
    def test_catvoxel_material_volumes_enabled(self, mock_rawwrite, mock_dumpjson, mock_printsetup, mock_makevol, mock_getmu, mock_feval):
        cfg = MagicMock()
        cfg.material_volumes = True
        cfg.write_vp = True
        cfg.make_img_kv = [100]
        cfg.spec.Evec = None
        cfg.phantom.callback = 'phantom_callback'
        cfg.phantom.Materials = ['Water', 'Bone']
        cfg.phantom.numberOfMaterials = 2
        cfg.phantom.projectorCallback = 'C_Projector_Analytic'
        cfg.phantom.filename = self.test_file
        cfg.recon.imageSize = 128
        cfg.recon.sliceCount = 64
        cfg.recon.fov = 256.0
        cfg.recon.sliceThickness = 1.0
        cfg.recon.centerOffset = [0.0, 0.0, 0.0]
        cfg.vol_os = 1
        dummy_volume = np.ones((2, 64, 128, 128), dtype=np.single)
        mock_makevol.return_value = dummy_volume
        catvoxel(cfg)
        mock_rawwrite.assert_called()
        mock_dumpjson.assert_called_once()
        mock_printsetup.assert_called_once()

    @patch('gecatsim.pyfiles.catvoxel.feval', return_value=MagicMock())
    @patch('gecatsim.pyfiles.catvoxel.GetMu', return_value=MagicMock())
    @patch('gecatsim.pyfiles.catvoxel.MakeAllVolumes')
    @patch('gecatsim.pyfiles.catvoxel.PrintPhantomSetup')
    @patch('gecatsim.pyfiles.catvoxel.dumpjson')
    @patch('builtins.open', new_callable=mock_open)
    def test_catvoxel_attenuation_volume(self, mock_file, mock_dumpjson, mock_printsetup, mock_makevol, mock_getmu, mock_feval):
        cfg = MagicMock()
        cfg.material_volumes = False
        cfg.write_vp = False
        cfg.make_img_kv = [120]
        cfg.spec.Evec = None
        cfg.phantom.callback = 'phantom_callback'
        cfg.phantom.Materials = ['Water']
        cfg.phantom.numberOfMaterials = 1
        cfg.phantom.projectorCallback = 'C_Projector_Analytic'
        cfg.phantom.filename = self.test_file
        cfg.recon.imageSize = 64
        cfg.recon.sliceCount = 32
        cfg.recon.fov = 128.0
        cfg.recon.sliceThickness = 1.0
        cfg.recon.centerOffset = [0.0, 0.0, 0.0]
        cfg.vol_os = 1
        dummy_volume = np.ones((1, 32, 64, 64), dtype=np.single)
        mock_makevol.return_value = dummy_volume
        catvoxel(cfg)
        mock_file.assert_called_once()

    def test_invalid_projector_callback(self):
        cfg = MagicMock()
        cfg.phantom = MagicMock()
        cfg.phantom.projectorCallback = 'Invalid_Callback'
        with self.assertRaises(ValueError):
            MakeAllVolumes(np.zeros((1,1,1,1), dtype=np.single), cfg, 1, False)

    @patch('gecatsim.pyfiles.catvoxel.feval', return_value=MagicMock())
    @patch('gecatsim.pyfiles.catvoxel.GetMu', return_value=MagicMock())
    @patch('gecatsim.pyfiles.catvoxel.MakeAllVolumes')
    @patch('gecatsim.pyfiles.catvoxel.PrintPhantomSetup')
    @patch('gecatsim.pyfiles.catvoxel.dumpjson')
    @patch('builtins.open', new_callable=mock_open)
    def test_missing_material_volumes_attribute(self, mock_file, mock_dumpjson, mock_printsetup, mock_makevol, mock_getmu, mock_feval):
        cfg = MagicMock()
        del cfg.material_volumes
        cfg.write_vp = False
        cfg.make_img_kv = [120]
        cfg.spec.Evec = None
        cfg.phantom.callback = 'phantom_callback'
        cfg.phantom.Materials = ['Water']
        cfg.phantom.numberOfMaterials = 1
        cfg.phantom.projectorCallback = 'C_Projector_Analytic'
        cfg.phantom.filename = self.test_file
        cfg.recon.imageSize = 64
        cfg.recon.sliceCount = 32
        cfg.recon.fov = 128.0
        cfg.recon.sliceThickness = 1.0
        cfg.recon.centerOffset = [0.0, 0.0, 0.0]
        cfg.vol_os = 1
        dummy_volume = np.ones((1, 32, 64, 64), dtype=np.single)
        mock_makevol.return_value = dummy_volume
        catvoxel(cfg)
        self.assertEqual(cfg.spec.Evec,  np.array(cfg.make_img_kv))

    @patch('gecatsim.pyfiles.catvoxel.feval', return_value=MagicMock())
    @patch('gecatsim.pyfiles.catvoxel.MakeAllVolumes')
    def test_zero_image_size(self, mock_makevol, mock_feval):
        cfg = MagicMock()
        cfg.recon.imageSize = 0
        cfg.recon.fov = 128.0
        cfg.recon.sliceThickness = 1.0
        cfg.recon.centerOffset = [0.0, 0.0, 0.0]
        cfg.recon.sliceCount = 32
        cfg.phantom.projectorCallback = 'C_Projector_Analytic'
        cfg.phantom.Materials = ['Water']
        cfg.phantom.numberOfMaterials = 1
        cfg.phantom.filename = self.test_file
        cfg.material_volumes = False
        cfg.write_vp = False
        cfg.make_img_kv = [120]
        cfg.spec.Evec = None
        cfg.phantom.callback = 'phantom_callback'
        cfg.vol_os = 1
        mock_feval.return_value = cfg
        with self.assertRaises(ZeroDivisionError):
            catvoxel(cfg)

    @patch('gecatsim.pyfiles.catvoxel.feval', return_value=MagicMock())
    @patch('gecatsim.pyfiles.catvoxel.GetMu', return_value=MagicMock())
    @patch('gecatsim.pyfiles.catvoxel.MakeAllVolumes')
    @patch('gecatsim.pyfiles.catvoxel.PrintPhantomSetup')
    @patch('builtins.open', new_callable=mock_open)
    def test_missing_spec_attribute(self, mock_file, mock_printsetup, mock_makevol, mock_getmu, mock_feval):
        cfg = MagicMock()
        del cfg.spec
        cfg.material_volumes = False
        cfg.make_img_kv = [120]
        cfg.phantom.callback = 'phantom_callback'
        cfg.phantom.Materials = ['Water']
        cfg.phantom.numberOfMaterials = 1
        cfg.phantom.projectorCallback = 'C_Projector_Analytic'
        cfg.phantom.filename = self.test_file
        cfg.recon.imageSize = 64
        cfg.recon.sliceCount = 32
        cfg.recon.fov = 128.0
        cfg.recon.sliceThickness = 1.0
        cfg.recon.centerOffset = [0.0, 0.0, 0.0]
        cfg.vol_os = 1
        dummy_volume = np.ones((1, 32, 64, 64), dtype=np.single)
        mock_makevol.return_value = dummy_volume
        with self.assertRaises(AttributeError):
            catvoxel(cfg)

    @patch('builtins.print')
    def test_print_phantom_setup_output(self, mock_print):
        cfg = MagicMock()
        cfg.Nx = 64
        cfg.Ny = 64
        cfg.Nz = 32
        cfg.dx = 1.0
        cfg.dy = 1.0
        cfg.dz = 1.0
        cfg.xoff = 0.0
        cfg.yoff = 0.0
        cfg.zoff = 0.0
        PrintPhantomSetup(cfg)
        self.assertEqual(mock_print.call_count, 9)

    def test_dumpjson_creates_valid_json(self):
        test_data = {"a": 1, "b": [2, 3]}
        test_path = os.path.join(self.test_dir.name, "test.json")
        dumpjson(test_path, test_data)
        with open(test_path, "r") as f:
            loaded = json.load(f)
        self.assertEqual(loaded, test_data)

    @patch('gecatsim.pyfiles.catvoxel.feval')
    def test_make_all_volumes_valid_callback(self, mock_feval):
        cfg = MagicMock()
        cfg.phantom.projectorCallback = 'C_Projector_Analytic'
        cfg.Nx, cfg.Ny, cfg.Nz = 64, 64, 32
        cfg.xoff, cfg.yoff, cfg.zoff = 0.0, 0.0, 0.0
        cfg.dx, cfg.dy, cfg.dz = 1.0, 1.0, 1.0
        cfg.vol_os = 1
        volume = np.zeros((1, 32, 64, 64), dtype=np.single)
        mock_feval.return_value = volume
        result = MakeAllVolumes(volume, cfg, 1, False)
        mock_feval.assert_called_once()
        np.testing.assert_array_equal(result, volume)

    @patch('gecatsim.pyfiles.catvoxel.feval', return_value=MagicMock())
    @patch('gecatsim.pyfiles.catvoxel.MakeAllVolumes')
    @patch('gecatsim.pyfiles.catvoxel.PrintPhantomSetup')
    @patch('gecatsim.pyfiles.catvoxel.xc.rawwrite')
    @patch('gecatsim.pyfiles.catvoxel.dumpjson')
    def test_missing_make_img_kv(self, mock_dumpjson, mock_rawwrite, mock_printsetup, mock_makevol, mock_feval):
        cfg = MagicMock()
        del cfg.make_img_kv
        cfg.material_volumes = False
        cfg.write_vp = False
        cfg.spec = MagicMock()
        cfg.phantom.callback = 'phantom_callback'
        cfg.phantom.Materials = ['Water']
        cfg.phantom.numberOfMaterials = 1
        cfg.phantom.projectorCallback = 'C_Projector_Analytic'
        cfg.phantom.filename = self.test_file
        cfg.recon.imageSize = 64
        cfg.recon.sliceCount = 32
        cfg.recon.fov = 128.0
        cfg.recon.sliceThickness = 1.0
        cfg.recon.centerOffset = [0.0, 0.0, 0.0]
        cfg.vol_os = 1
        dummy_volume = np.ones((1, 32, 64, 64), dtype=np.single)
        mock_makevol.return_value = dummy_volume
        catvoxel(cfg)
        self.assertTrue((cfg.spec.Evec == np.array([120])).all())

    @patch('gecatsim.pyfiles.catvoxel.feval', return_value=MagicMock())
    @patch('gecatsim.pyfiles.catvoxel.MakeAllVolumes')
    @patch('gecatsim.pyfiles.catvoxel.PrintPhantomSetup')
    @patch('gecatsim.pyfiles.catvoxel.xc.rawwrite')
    @patch('gecatsim.pyfiles.catvoxel.dumpjson')
    def test_material_volumes_without_write_vp(self, mock_dumpjson, mock_rawwrite, mock_printsetup, mock_makevol, mock_feval):
        cfg = MagicMock()
        cfg.material_volumes = True
        cfg.write_vp = False
        cfg.make_img_kv = [100]
        cfg.spec = MagicMock()
        cfg.phantom.callback = 'phantom_callback'
        cfg.phantom.Materials = ['Water', 'Bone']
        cfg.phantom.numberOfMaterials = 2
        cfg.phantom.projectorCallback = 'C_Projector_Analytic'
        cfg.phantom.filename = self.test_file
        cfg.recon.imageSize = 64
        cfg.recon.sliceCount = 32
        cfg.recon.fov = 128.0
        cfg.recon.sliceThickness = 1.0
        cfg.recon.centerOffset = [0.0, 0.0, 0.0]
        cfg.vol_os = 1
        dummy_volume = np.ones((2, 32, 64, 64), dtype=np.single)
        mock_makevol.return_value = dummy_volume

        # Mock feval to return cfg unchanged
        mock_feval.side_effect = lambda cb, c: c

        catvoxel(cfg)
        mock_rawwrite.assert_not_called()
        mock_dumpjson.assert_not_called()

if __name__ == '__main__':
    unittest.main()