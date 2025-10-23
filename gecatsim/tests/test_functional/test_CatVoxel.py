import json
import os
import unittest.mock

import gecatsim as xc
from gecatsim.pyfiles.catvoxel import catvoxel


class Test_CatVoxel(unittest.TestCase):

    def initialize(self):
        ct = xc.CatSim()
        cfg = ct.get_current_cfg()
        # -------------------
        # some common settings

        # Make separate volumes of densities (in gm/cm^3) for each materia
        cfg.material_volumes = 1
        # oversampling in making volumes
        cfg.vol_os = 11
        # Sets the kV when cfg.material_volumes is set to 0
        # in this case, the attenuation coefficients at this kVp will be saved
        # cfg.make_img_kv = 120
        # save results
        cfg.write_vp = 1

        cfg.recon.imageSize = 512
        cfg.recon.sliceCount = 10
        cfg.recon.sliceThickness = 0.568
        cfg.phantom.scale = 1.0
        return cfg

    def test_catvoxel_Analytic_Phantom(self):
        cfg = self.initialize()

        self.base_name = "CTDI_16cm_WaterAirPEBoneChambers"
        self.expected_files = [
            f"{self.base_name}.json",
            f"{self.base_name}.VolumeFraction_pmma",
            f"{self.base_name}.VolumeFraction_water",
            f"{self.base_name}.VolumeFraction_air",
            f"{self.base_name}.VolumeFraction_polyethylene",
            f"{self.base_name}.VolumeFraction_ICRU_skeleton_cortical_bone_adult"
        ]

        # analytic phantom
        print("running catvoxel for analytic phantom...")
        cfg.phantom.projectorCallback = 'C_Projector_Analytic'
        cfg.phantom.callback = 'Phantom_Analytic'
        cfg.phantom.filename = 'CTDI_16cm_WaterAirPEBoneChambers.ppm'
        cfg.recon.fov = 500.0
        cfg.recon.centerOffset = [0, 0, 0]
        catvoxel(cfg)

        self.test_output_files_exist()

        data = {}

        with open(f"{self.base_name}.json", "r") as f:
            data = json.load(f)

        ##test keys
        expected_keys = [
            "n_materials", "mat_name", "volumefractionmap_filename",
            "volumefractionmap_datatype", "cols", "rows", "slices",
            "x_size", "y_size", "z_size", "x_offset", "y_offset", "z_offset"
        ]
        for key in expected_keys:
            self.assertIn(key, data, f"Missing key: {key}")

        ##test material names
        expected_materials = [
            "pmma", "water", "air", "polyethylene", "ICRU_skeleton_cortical_bone_adult"
        ]
        self.assertEqual(data["mat_name"], expected_materials)

        ## test volume file names
        for mat in data["mat_name"]:
            expected_filename = f"CTDI_16cm_WaterAirPEBoneChambers.VolumeFraction_{mat}"
            self.assertIn(expected_filename, data["volumefractionmap_filename"])

        ## test dimension and offset
        self.assertTrue(all(v == 512 for v in data["cols"]), "Unexpected cols")
        self.assertTrue(all(v == 512 for v in data["rows"]), "Unexpected rows")
        self.assertTrue(all(v == 10 for v in data["slices"]), "Unexpected slices")
        self.assertTrue(all(v == 0.9765625 for v in data["x_size"]), "Unexpected x_size")
        self.assertTrue(all(v == 0.9765625 for v in data["y_size"]), "Unexpected y_size")
        self.assertTrue(all(v == 0.568 for v in data["z_size"]), "Unexpected z_size")
        self.assertTrue(all(v == 256.5 for v in data["x_offset"]), "Unexpected x_offset")
        self.assertTrue(all(v == 256.5 for v in data["y_offset"]), "Unexpected y_offset")
        self.assertTrue(all(v == 5.5 for v in data["z_offset"]), "Unexpected z_offset")


    def test_catvoxel_NCAT_Phantom(self):
        cfg = self.initialize()

        self.base_name = "vmale50_chest_less_surfaces"
        self.expected_files = [
            f"{self.base_name}.json",
            f"{self.base_name}.VolumeFraction_ncat_adipose"
        ]

        # NCAT phantom
        print("running catvoxel for NCAT phantom...")
        cfg.phantom.projectorCallback = 'C_Projector_NCAT'
        cfg.phantom.callback = 'Phantom_NCAT'
        cfg.phantom.filename = 'vmale50_chest_less_surfaces.nrb'
        cfg.recon.fov = 400.0
        cfg.recon.centerOffset = [0, 0, 0]
        catvoxel(cfg)

        self.test_output_files_exist()

        # Validate JSON contents
        with open(f"{self.base_name}.json", "r") as f:
            data = json.load(f)

        self.assertEqual(data["n_materials"], 1)
        self.assertEqual(data["mat_name"], ["ncat_adipose"])
        self.assertEqual(data["volumefractionmap_filename"], [f"{self.base_name}.VolumeFraction_ncat_adipose"])
        self.assertEqual(data["volumefractionmap_datatype"], ["float"])
        self.assertEqual(data["cols"], [512])
        self.assertEqual(data["rows"], [512])
        self.assertEqual(data["slices"], [10])
        self.assertEqual(data["x_size"], [0.78125])
        self.assertEqual(data["y_size"], [0.78125])
        self.assertEqual(data["z_size"], [0.568])
        self.assertEqual(data["x_offset"], [256.5])
        self.assertEqual(data["y_offset"], [256.5])
        self.assertEqual(data["z_offset"], [5.5])

    def test_catvoxel_Polygonal_Phantom(self):
        cfg = self.initialize()

        self.base_name = "female_adult_average_lung_lesions_reduced"
        self.expected_files = [
            f"{self.base_name}.json",
            f"{self.base_name}.VolumeFraction_ncat_muscle"
        ]

        # polygonal phantom
        print("running catvoxel for polygonal phantom...")
        cfg.phantom.projectorCallback = 'C_Projector_Polygon'
        cfg.phantom.callback = 'Phantom_Polygonal'
        cfg.phantom.filename = 'female_adult_average_lung_lesions_reduced.nrb'
        cfg.recon.fov = 50.0
        cfg.recon.centerOffset = [-57, -51, 102]
        catvoxel(cfg)

        self.test_output_files_exist()

        # Check output files exist
        for fname in self.expected_files:
            self.assertTrue(os.path.exists(fname), f"Missing output file: {fname}")

        # Validate JSON contents
        with open(f"{self.base_name}.json", "r") as f:
            data = json.load(f)

        self.assertEqual(data["n_materials"], 1)
        self.assertEqual(data["mat_name"], ["ncat_muscle"])
        self.assertEqual(data["volumefractionmap_filename"], [f"{self.base_name}.VolumeFraction_ncat_muscle"])
        self.assertEqual(data["volumefractionmap_datatype"], ["float"])
        self.assertEqual(data["cols"], [512])
        self.assertEqual(data["rows"], [512])
        self.assertEqual(data["slices"], [10])
        self.assertEqual(data["x_size"], [0.09765625])
        self.assertEqual(data["y_size"], [0.09765625])
        self.assertEqual(data["z_size"], [0.568])
        self.assertEqual(data["x_offset"], [840.18])
        self.assertEqual(data["y_offset"], [778.74])
        self.assertEqual(data["z_offset"], [-174.0774647887324])

    def test_output_files_exist(self):
        for fname in self.expected_files:
            self.assertTrue(os.path.exists(fname), f"Missing output file: {fname}")
