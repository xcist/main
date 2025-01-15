import unittest
from unittest.mock import patch
import numpy as np
import os
from gecatsim.pyfiles.writeHdf5 import writeHdf5, getNumExtraViews, calcImgTableParams
from gecatsim.pyfiles.CommonTools import CFG

class TestWriteHdf5(unittest.TestCase):
    def setUp(self):
        self.cfg = CFG("../examples/cfg/Phantom_Sample",
                       "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")
        # Add missing attributes
        self.cfg.recon_slice_thickness = 1.0
        self.cfg.recon_mu = 0.02
        self.cfg.col_size = 1.0
        self.cfg.sid = 1000.0
        self.cfg.sdd = 1500.0
        self.cfg.results_basename = "test"
        self.cfg.total_n_viewsverbo = 10
        self.cfg.col_count = 512
        self.cfg.row_count = 512
        self.cfg.betas = np.linspace(0, 2*np.pi, 10)
        self.cfg.mA = 200
        self.cfg.rotation_period = 1.0
        self.cfg.views_per_rotation = 360
        self.cfg.gammas = np.array([[0, 1]])
        self.cfg.sbd = 1.0
        self.cfg.col_offset = 0
        self.cfg.start_view = 0
        self.cfg.total_n_views = 360
        self.cfg.start_z = 0.0
        self.cfg.recon_zcenter = 0.0
        self.cfg.recon_planes = 10
        self.cfg.table_speed = 1.0

    @patch("h5py.File", autospec=True)
    def test_writeHdf5(self, mock_h5py):
        # Create a temporary file to use as the input file
        with open("test.prep", "wb") as f:
            f.write(b'\x00' * 512 * 512 * 4 * self.cfg.total_n_viewsverbo)

        # Verify the file content
        with open("test.prep", "rb") as f:
            data = f.read()
            print(f"File size: {len(data)} bytes")

        # Create the output folder if it doesn't exist
        if not os.path.exists("output_folder"):
            os.makedirs("output_folder")

        writeHdf5(self.cfg, "output_folder")
        mock_h5py.assert_called_once_with("test.h5", "w")

        # Clean up the temporary file
        os.remove("test.prep")

    def test_getNumExtraViews(self):
        nExtraViews = getNumExtraViews(self.cfg)
        self.assertEqual(nExtraViews, 12)

    def test_calcImgTableParams(self):
        imgTable = calcImgTableParams(self.cfg, 1.0)
        self.assertEqual(imgTable["number_of_images"], 0)

if __name__ == '__main__':
    unittest.main()