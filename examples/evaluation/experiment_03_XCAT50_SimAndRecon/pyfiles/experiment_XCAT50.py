# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

# Author: Paul FitzGerald
# Date: April 27, 2022
#
# Purpose: This is an XCIST "experiment file" that is used to evaluate several aspects of XCIST simulation and recon
# using an XCAT phantom. The default config files are used for everything except the phantom - for that, you will need
# Phantom_XCAT50.cfg, which should be included in the same folder where you found this file.
#
# This "experiment file" performs several experiments, and each experiment requires several simulations and
# reconstructions. Results are written to the folder defined by the system environment variable XCIST_UserPath.
# Sub-folders are for each sim/recon. DON'T make this path point to a location within your local copy of the CatSim
# repository, or the output files (which are numerous and large) might end up in the repository!
#
# The experiments performed evaluate image quality (IQ) related to the apparent spatial resolution
# as two relevant parametes are varied: focal spot (FS) size and detector column width. The experiments are:
# 01_xx through 03_xx. focal spot (FS) size
# xx_01 through xx_03 Vary detector column width
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

    string1 = "focalspotWidth = {}, srcXSampleCount = {};".format(                                                              \
              cfg.scanner.focalspotWidth, cfg.physics.srcXSampleCount)
    string2 = "detectorColSize = {}, detectorColCount = {}, colSampleCount = {};".format(                                       \
              cfg.scanner.detectorColSize, cfg.scanner.detectorColCount, cfg.physics.colSampleCount)
    string3 = "viewsPerRotation = {}, viewCount = {}, stopViewId = {}, viewSampleCount = {}.".format(                           \
              cfg.protocol.viewsPerRotation, cfg.protocol.viewCount, cfg.protocol.stopViewId, cfg.physics.viewSampleCount)
    return WLString(cfg) + "\n" + string1 + "\n" + string2 + "\n" + string3


##--------- Initialize

userPath = getUserPath()
xc.CommonTools.my_path.add_search_path(userPath)

# Use the default cfg parameters found in the default .cfg files, except use a specific phantom file.

cfg = xc.CatSim(xc.CommonTools.my_path.find("cfg", "Phantom_XCAT50.cfg", ""))
cfg.experimentDirectory = os.path.join(userPath, "examples", "evaluation", "experiment_03_XCAT50_SimAndRecon")

# These are changes to the default config parameters to be used for this experiment.

cfg.scanner.detectorRowsPerMod = 1
cfg.scanner.detectorRowCount = cfg.scanner.detectorRowsPerMod
cfg.scanner.detectorColOffset = 0.25
# cfg.scanner.focalspotWidth                # Will be revised in the loop below.
# cfg.scanner.detectorColSize               # Will be revised in the loop below.
# cfg.scanner.detectorColCount              # Will be revised in the loop below.

cfg.protocol.mA = 1200
cfg.protocol.spectrumScaling = 0.931
cfg.protocol.bowtie = "small.txt"
# cfg.protocol.viewsPerRotation             # Will be revised in the loop below.
# cfg.protocol.viewCount                    # Will be revised in the loop below.
# cfg.protocol.stopViewId                   # Will be revised in the loop below.

cfg.physics.srcZSampleCount = 1
cfg.physics.rowSampleCount = 1
# cfg.physics.srcXSampleCount               # Will be revised in the loop below.
# cfg.physics.colSampleCount                # Will be revised in the loop below.
# cfg.physics.viewSampleCount               # Will be revised in the loop below.

cfg.recon.fov = 60.0
cfg.recon.centerOffset = [0.0, 50.0, 0.0]
cfg.recon.kernelType = 'Bone'
cfg.recon.displayWindowMin = -400       # In HU.
cfg.recon.displayWindowMax = 400        # In HU.
cfg.recon.displayWindow = cfg.recon.displayWindowMax - cfg.recon.displayWindowMin
cfg.recon.displayLevel = (cfg.recon.displayWindowMax + cfg.recon.displayWindowMin)/2
# For publication.
cfg.recon.saveImagePictureFiles = True
# For test development - comment out is desired for published figures.
# cfg.recon.displayImagePictureAxes = False
# cfg.recon.displayImagePictureTitles = False
# For single-row simulations, use the native slice thickness.
if cfg.scanner.detectorRowCount == 1:
    cfg.recon.sliceThickness = cfg.scanner.detectorRowSize*cfg.scanner.sid/cfg.scanner.sdd

# Top-level cfg parameters related to control of this experiment.
cfg.waitForKeypress = False             # Wait for keypress after plot display?
cfg.do_Sim = False                      # The simulation is usually run except when only varying recon parameters.
cfg.displaySimProjectionPlots = False   # Flag to display plots to screen.
cfg.do_Recon = True                     # The recon is usually run except when only varying display parameters.

##--------- Loop through experiment parameters

focalspotWidth_settings   = [[4.0,  4.0,  4.0 ], [0.25, 0.25, 0.25]]
srcXSampleCount_settings  = [[ 17,   17,   17 ], [ 3  ,  3  ,  3  ]]
detectorColSize_settings  = [[2.0 , 1.0 , 0.25], [2.0 , 1.0 , 0.25]]
colSampleCount_settings   = [[11  ,  7  ,  3  ], [11  ,  7  ,  3  ]]
viewsPerRotation_settings = [[1000, 1000, 4000], [1000, 1000, 4000]]
viewSampleCount_settings  = [[ 2  ,  2  ,  1  ], [ 2  ,  2  ,  1  ]]

for focalspotWidth_Index in range(0, len(focalspotWidth_settings)):
    for detectorColSize_Index in range(0, len(detectorColSize_settings[0])):
        cfg.scanner.focalspotWidth = focalspotWidth_settings[focalspotWidth_Index][detectorColSize_Index]
        cfg.physics.srcXSampleCount = srcXSampleCount_settings[focalspotWidth_Index][detectorColSize_Index]
        cfg.scanner.detectorColSize = detectorColSize_settings[focalspotWidth_Index][detectorColSize_Index]
        cfg.scanner.detectorColCount = int(480/cfg.scanner.detectorColSize)
        cfg.physics.colSampleCount = colSampleCount_settings[focalspotWidth_Index][detectorColSize_Index]
        cfg.protocol.viewsPerRotation = viewsPerRotation_settings[focalspotWidth_Index][detectorColSize_Index]
        cfg.protocol.viewCount = cfg.protocol.viewsPerRotation
        cfg.protocol.stopViewId = cfg.protocol.viewCount - 1
        cfg.physics.viewSampleCount = viewSampleCount_settings[focalspotWidth_Index][detectorColSize_Index]
        cfg.experimentName = "{:02d}_{:02d}_focalspotWidth{:4.2f}_detectorColSize{:4.2f}".format(     \
                             focalspotWidth_Index + 1, detectorColSize_Index + 1,                     \
                             cfg.scanner.focalspotWidth, cfg.scanner.detectorColSize)
        cfg.experimentName = cfg.experimentName.replace(".", "p")

        # Initialize the experiment directory.
        cfg = initializeExperimentDirectory(cfg)

        # Get the title for the recon images.
        cfg.reconImageTitle = getReconImageTitle(cfg)

        # Do the experiment.
        doExperiment(cfg)
