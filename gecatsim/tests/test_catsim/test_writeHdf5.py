import unittest
from unittest.mock import patch, MagicMock
import numpy as np
import os
from gecatsim.pyfiles.writeHdf5 import writeHdf5, getNumExtraViews, calcImgTableParams
from gecatsim.pyfiles.CommonTools import CFG

cfg = CFG("../examples/cfg/Phantom_Sample",
          "../examples/cfg/Scanner_Sample_generic",
          "../examples/cfg/Protocol_Sample_axial")

cfg.recon_slice_thickness = 1.0
cfg.recon_mu = 0.02
cfg.col_size = 1.0
cfg.sid = 1000.0
cfg.sdd = 1500.0
cfg.results_basename = "test"
cfg.total_n_viewsverbo = 10
cfg.col_count = 512
cfg.row_count = 512
cfg.betas = np.linspace(0, 2 * np.pi, 10)
cfg.mA = 200
cfg.rotation_period = 1.0
cfg.views_per_rotation = 360
cfg.gammas = np.array([[0, 1]])
cfg.sbd = 1.0
cfg.col_offset = 0
cfg.start_view = 0
cfg.total_n_views = 360
cfg.start_z = 0.0
cfg.recon_zcenter = 0.0
cfg.recon_planes = 10
cfg.table_speed = 1.0

@patch("h5py.File", autospec=True)
def test_writeHdf5(mock_h5py):
    # Create a temporary file to use as the input file
    with open("test.prep", "wb") as f:
        f.write(b'\x00' * 512 * 512 * 4 * cfg.total_n_viewsverbo)

    # Verify the file content
    with open("test.prep", "rb") as f:
        data = f.read()
        print(f"File size: {len(data)} bytes")

    # Create the output folder if it doesn't exist
    if not os.path.exists("output_folder"):
        os.makedirs("output_folder")

    writeHdf5(cfg, "output_folder")
    mock_h5py.assert_called_once_with("test.h5", "w")

    # Clean up the temporary file
    os.remove("test.prep")

def test_getNumExtraViews():
    nExtraViews = getNumExtraViews(cfg)
    assert nExtraViews == 12

def test_calcImgTableParams():
    imgTable = calcImgTableParams(cfg, 1.0)
    assert imgTable["number_of_images"] == 0

if __name__ == '__main__':
    unittest.main()