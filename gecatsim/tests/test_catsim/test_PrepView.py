import os
import unittest.mock

import gecatsim as xc


class Test_PrepView(unittest.TestCase):
    def test_prepView(self):

        ct = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")  # initialization

        file_names = ["test.air", "test.offset", "test.prep", "test.scan"]

        ##--------- Make changes to parameters (optional)
        ct.resultsName = "test"
        ct.protocol.viewsPerRotation = 50
        ct.protocol.viewCount = ct.protocol.viewsPerRotation
        ct.protocol.stopViewId = ct.protocol.viewCount - 1

        ##--------- Run simulation
        ct.run_all()

        assert (os.path.exists('test.prep') == True)

        for file_name in file_names:
            if os.path.exists(file_name):
                os.remove(file_name)
