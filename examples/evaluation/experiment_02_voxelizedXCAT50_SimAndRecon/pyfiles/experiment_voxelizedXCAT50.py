# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

# Author: Paul FitzGerald
# Date: April 18, 2022
#
# Purpose: This is an XCIST "experiment file" that is used to evaluate several aspects of XCIST simulation and recon
# using a voxelized version of an XCAT phantom. The default config files are used for everything except the phantom - for
# that, you will need Phantom_voxelizedXCAT50.cfg, which should be included in the same folder where you found this file.
#
# This "experiment file" performs several experiments, and each experiment requires several simulations and
# reconstructions. Results are written to the folder defined by the system environment variable XCIST_UserPath.
# Sub-folders are for each sim/recon. DON'T make this path point to a location within your local copy of the CatSim
# repository, or the output files (which are numerous and large) might end up in the repository!
#
# The experiments performed evaluate image quality (IQ) related to the view aliasing and noise as two relevant
# parameters are varied: views per rotation (V/R) and radiation dose via tube current (mA). The experiments are:
# 01_xx through 03_xx. Vary tube current (mA)
# xx_01 through xx_04 Vary views/rotation
#
# The overall process is:
# a. Define "base" config (see "##--------- Initialize").
# b. Define parameters to run (see "##--------- Loop through experiment parameters").
# c. Loop through all the sim/recons defined in b, and for each sim/recon:
# d.   Define the specific parameters for the sim/recon.
# e.   Get the title for the recon images - function getReconImageTitle().
#      This includes parameters relevant to the experiment.
# f.   Run the simulation - function runSim().
# g.   Run the recon - function runRecon().
#
# Note that if you set
#   cfg.do_Sim = False
# the simulation will not run, but there needs to be a valid simulation result based on the defined simulation parameters,
# and the projection plots will be desplayed from that simulation.
# Also, if you set
#   cfg.do_Recon = False
# the reconstruction will not run, but there needs to be a valid reconstruction result based on the defined simulation and
# reconstuction parameters, and the reconstruction results will be displayed/saved based on the reconstruction display
# parameters such as window/level and titles.

import os
import catsim as xc
from my_commonTools import *


def getReconImageTitle(cfg):

    # Most experiments have only one image, so don't add that to the title.
    cfg.addSliceInfoToReconImageTitle = False

    formatString = "cfg.protocol.mA = {}; cfg.protocol.ciewsPerRotation = {};"
    string1 = formatString.format(cfg.protocol.mA, cfg.protocol.viewsPerRotation)
    return WLString(cfg) + "\n" + string1


##--------- Initialize

userPath = getUserPath()
xc.CommonTools.my_path.add_search_path(userPath)

# Use the default cfg parameters found in the default .cfg files, except use a specific phantom file.

cfg = xc.CatSim(xc.CommonTools.my_path.find("cfg", "Phantom_voxelizedXCAT50.cfg", ""))
cfg.phantomGender = "female"
cfg.experimentDirectory = os.path.join(userPath, "examples", "evaluation", "experiment_02_voxelizedXCAT50_SimAndRecon", cfg.phantomGender)
cfg.phantom.filename = "adult_" + cfg.phantomGender + "_50percentile_chest_slab_400.json"

# These are changes to the default config parameters to be used for this experiment.

low_mA = False

cfg.scanner.detectorRowsPerMod = 1
cfg.scanner.detectorRowCount = cfg.scanner.detectorRowsPerMod
cfg.scanner.detectorColOffset = 0.25
if low_mA:
    pass                            # if the mA is super-low, use the default eNoise
else:
    cfg.scanner.eNoise = 20000.0    # If the mA is high enough, need to boost eNoise to see the effect.

cfg.protocol.spectrumScaling = 0.931
# cfg.protocol.mA                           # Will be revised in the loop below.
# cfg.protocol.viewsPerRotation             # Will be revised in the loop below.
# cfg.protocol.viewCount                    # Will be revised in the loop below.
# cfg.protocol.stopViewId                   # Will be revised in the loop below.
cfg.protocol.maxPrep = -1       # To be sure that streak artifacts due to small values from e-noise are not supressed.

cfg.physics.energyCount = 12 # Using 120 kVp, so 10 kV/bin.
cfg.physics.srcXSampleCount = 2
cfg.physics.srcZSampleCount = 1
cfg.physics.rowSampleCount = 1
cfg.physics.colSampleCount = 1
cfg.physics.viewSampleCount = 1

cfg.recon.fov = 325.0                   # in mm.
cfg.recon.displayWindowMin = -350       # In HU.
cfg.recon.displayWindowMax = 50        # In HU.
cfg.recon.displayWindow = cfg.recon.displayWindowMax - cfg.recon.displayWindowMin
cfg.recon.displayLevel = (cfg.recon.displayWindowMax + cfg.recon.displayWindowMin)/2
# For publication.
cfg.recon.saveImagePictureFiles = True
# For test development - comment out if desired for published figures.
# cfg.recon.displayImagePictureAxes = False
# cfg.recon.displayImagePictureTitles = False
# For single-row simulations, use the native slice thickness.
if cfg.scanner.detectorRowCount == 1:
    cfg.recon.sliceThickness = cfg.scanner.detectorRowSize*cfg.scanner.sid/cfg.scanner.sdd

# Top-level cfg parameters related to control of this experiment.
cfg.waitForKeypress = False             # Wait for keypress after plot display?
cfg.do_Sim = False                      # The simulation is usually run except when only varying recon parameters.
cfg.displaySimProjectionPlots = False   # Flag to display plots to screen.
cfg.do_Recon = False                    # The recon is usually run except when only varying display parameters.

##--------- Loop through experiment parameters

viewsPerRotation_settings = [128, 256, 512, 1024]
if low_mA:
    mA_settings = [2, 10, 50]
else:
    mA_settings = [16, 100, 500]

for mA_Index in range(0, len(mA_settings)):
    for viewsPerRotation_Index in range(0, len(viewsPerRotation_settings)):
        cfg.protocol.mA = mA_settings[mA_Index]
        cfg.protocol.viewsPerRotation = viewsPerRotation_settings[viewsPerRotation_Index]
        cfg.protocol.viewCount = cfg.protocol.viewsPerRotation
        cfg.protocol.stopViewId = cfg.protocol.viewCount - 1
        cfg.experimentName = "{:02d}_{:02d}_{}mA_{}viewsPerRotation".format(     \
                             mA_Index + 1, viewsPerRotation_Index + 1,           \
                             cfg.protocol.mA, cfg.protocol.viewsPerRotation)

        # Initialize the experiment directory.
        cfg = initializeExperimentDirectory(cfg)

        # Get the title for the recon images.
        cfg.reconImageTitle = getReconImageTitle(cfg)

        # Do the experiment.
        doExperiment(cfg)
