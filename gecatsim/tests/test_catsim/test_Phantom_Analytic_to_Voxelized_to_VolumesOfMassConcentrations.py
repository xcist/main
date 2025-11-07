import unittest
import numpy as np
from unittest.mock import patch
import gecatsim as xc
from gecatsim.pyfiles.Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations import (
    Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations
)

class TestPhantomAnalyticToVoxelizedToVolumesOfMassConcentrations(unittest.TestCase):

    @patch('gecatsim.pyfiles.Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.catvoxel')
    @patch('os.path.getmtime')
    @patch('os.path.exists')
    def test_Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations(self, mock_exists, mock_getmtime, mock_catvoxel):

        ct = xc.CatSim()
        cfg = ct.get_current_cfg()

        # Set phantom parameters
        cfg.phantom.filename = 'CTDI_16cm_WaterAirPEBoneChambers.ppm'
        cfg.phantom.samples_xy = 128
        cfg.phantom.samples_voxelsize = 1.0
        cfg.phantom.samples_z = 64
        cfg.force_phantom_conversion = True
        cfg.write_vp = 1
        cfg.material_volumes = 1

        # Patch file system
        mock_exists.return_value = False
        mock_getmtime.return_value = 0

        mock_catvoxel.return_value = (
            np.zeros((2, 64, 128, 128), dtype=np.float32),
            ['water', 'bone']
        )

        Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations(cfg)

        self.assertEqual(cfg.phantom_filename, 'CTDI_16cm_WaterAirPEBoneChambers.vp')
        mock_catvoxel.assert_called_once()

if __name__ == '__main__':
    unittest.main()
