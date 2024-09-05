import os
import unittest.mock
from unittest.mock import patch

import numpy
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.WriteRawViewChunk import WriteRawViewChunk
import gecatsim as xc
from io import BytesIO, StringIO


class TestWriteRawViewChunk(unittest.TestCase):

    def test_write_raw_view_chunk_writes_file_in_writemode(self):
        ct = xc.CatSim()  # initialization

        ##--------- Make changes to parameters (optional)
        ct.resultsName = "test_sc"

        ct.sim = emptyCFG
        ct.sim.startViewId = 0
        ct.sim.isAirScan = True
        ct.sim.isOffsetScan = True
        ct.sim.isPhantomScan  = True

        ct.scanner.detectorRowsPerMod = 64
        ct.scanner.detectorRowCount = ct.scanner.detectorRowsPerMod
        ct.thisView = np.random.random([900, 18, 64])

        WriteRawViewChunk(ct, 0)

        scanTypeInd = [ct.sim.isAirScan, ct.sim.isOffsetScan, ct.sim.isPhantomScan].index(1)
        extName = ['.air', '.offset', '.scan'][scanTypeInd]
        fname = ct.resultsName + extName

        assert os.access(fname, os.W_OK)

        if os.path.exists(fname):
            os.remove(fname)

    def test_write_raw_view_chunk_writes_file_in_append(self):
        ct = xc.CatSim()  # initialization

        ##--------- Make changes to parameters (optional)
        ct.resultsName = "test_sc"

        ct.sim = emptyCFG
        ct.sim.startViewId = 0
        ct.sim.isAirScan = True
        ct.sim.isOffsetScan = True
        ct.sim.isPhantomScan = True
        ct.sim.stopViewId = 1

        ct.scanner.detectorRowsPerMod = 64
        ct.scanner.detectorRowCount = ct.scanner.detectorRowsPerMod
        ct.thisView = np.random.random([900, 18, 64])
        ct.dump_views = b''

        WriteRawViewChunk(ct, 1)

        scanTypeInd = [ct.sim.isAirScan, ct.sim.isOffsetScan, ct.sim.isPhantomScan].index(1)
        extName = ['.air', '.offset', '.scan'][scanTypeInd]
        fname = ct.resultsName + extName

        assert os.access(fname, os.R_OK)
        assert os.access(fname, os.W_OK)

        if os.path.exists(fname):
            os.remove(fname)