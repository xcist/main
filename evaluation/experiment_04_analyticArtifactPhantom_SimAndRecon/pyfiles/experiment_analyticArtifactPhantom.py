# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

# Author: Paul FitzGerald
# Date: May 16, 2022
#
# Purpose: This is an XCIST "experiment file" that is used to evaluate several aspects of XCIST simulation and recon
# using the artifacts phantom. The default config files are used for everything except the phantom - for that, you
# will need Phantom_analyticArtifactsPhantom.cfg, which should be included in the same folder where you found this file.
#
# The artifacts logo phantom used in this experiment is an analytic phantom. Therefore, the simulation portion of this
# evaluation only evaluates the analytic projector. However, the recon aspects of this experiment are independent of
# the phantom/projector type.
#
# This "experiment file" performs several experiments, and each experiment requires several simulations and
# reconstructions. Results are written to the folder defined by the system environment variable XCIST_UserPath.
# Sub-folders are for each sim/recon. DON'T make this path point to a location within your local copy of the CatSim
# repository, or the output files (which are numerous and large) might end up in the repository!
#
# The experiments performed include:
# 01. Artifacts simulation evaluation (6 simulations/reconstructions)
#
# Each sim/recon is independent.
# Each sim/recon is included in a list of sim/recons to run - see "##--------- Define experiment names".
# Each can be run or not by uncommenting or commenting them.
#
# The overall process is:
# a. Define "base" config (see "##--------- Initialize").
# b. Define a list of sim/recons to run (see "##--------- Define experiment names").
# c. Loop through all the sim/recons defined in b, and for each sim/recon:
# d.   Define the specific parameters for the sim/recon - function setExperimentParameters().
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
import copy
import catsim.pyfiles as catsim
from my_commonTools import *

def setExperimentParameters(cfg):

    experimentName = cfg.experimentName

    # Some experiments use fewer views.
    if experimentName == "01_02_ViewAliasing"           \
    or experimentName == "01_06_AllArtifacts":
        cfg.protocol.viewsPerRotation = 360
        cfg.protocol.viewCount = cfg.protocol.viewsPerRotation
        cfg.protocol.stopViewId = cfg.protocol.viewCount - 1
    
    # Some experiments use a polyenergetic spectrum and its estimated water mu.
    if experimentName == "01_03_BeamHardening"          \
    or experimentName == "01_06_AllArtifacts":
        cfg.physics.monochromatic = -1
        cfg.recon.mu = 0.02

    # Some experiments use scatter, but we need a 64-slice simulation for the scatter model to be correct.
    if experimentName == "01_04_Scatter"                \
    or experimentName == "01_06_AllArtifacts":
        cfg.scanner.detectorRowsPerMod = 64
        cfg.scanner.detectorRowCount = cfg.scanner.detectorRowsPerMod
        cfg.physics.scatterCallback = "Scatter_ConvolutionModel"
        cfg.physics.scatterKernelCallback = ""
        cfg.physics.scatterScaleFactor = 1
    
    # Some experiments use electonic and quantum noise.
    if experimentName == "01_05_Noise"                  \
    or experimentName == "01_06_AllArtifacts":
        cfg.physics.enableQuantumNoise = 1
        cfg.physics.enableElectronicNoise = 1

    return cfg


def getReconImageTitle(cfg):

    # Most experiments have only one image, so don't add that to the title.
    cfg.addSliceInfoToReconImageTitle = False

    return ""


##--------- Initialize

userPath = getUserPath()
catsim.CommonTools.my_path.add_search_path(userPath)

# Use the default cfg parameters found in the default .cfg files, except use a specific phantom file.

phantomCfgPathname = catsim.CommonTools.my_path.find("cfg", "Phantom_analyticArtifactPhantom.cfg", "")
cfg = catsim.CatSim.CatSim(phantomCfgPathname)
cfg.experimentDirectory = os.path.join(userPath, "examples", "evaluation", "experiment_04_analyticArtifactPhantom_SimAndRecon")

# These are changes to the defaults config parameters to be used for the "base" experiment,
# but some of these and some other defaults are overridden for certain experiments.

cfg.scanner.detectorRowsPerMod = 1                              # Will be revised for certain experiments.
cfg.scanner.detectorRowCount = cfg.scanner.detectorRowsPerMod   # Will be revised for certain experiments.
cfg.scanner.detectorColOffset = 0.25                            # Will be revised for certain experiments.

cfg.protocol.viewsPerRotation = 2000                            # Will be revised for certain experiments.
cfg.protocol.viewCount = cfg.protocol.viewsPerRotation          # Will be revised for certain experiments.
cfg.protocol.stopViewId = cfg.protocol.viewCount - 1            # Will be revised for certain experiments.
cfg.protocol.spectrumScaling = 0.931
cfg.protocol.mA = 500
cfg.protocol.maxPrep = -1

cfg.physics.energyCount = 12                                    # Using 120 kVp, so 10 kV/bin.
cfg.physics.monochromatic = 88                                  # Will be revised for certain experiments.
cfg.physics.colSampleCount = 4
cfg.physics.rowSampleCount = 1
cfg.physics.srcXSampleCount = 4
cfg.physics.srcYSampleCount = 1
cfg.physics.viewSampleCount = 1
cfg.physics.enableQuantumNoise = 0                              # Will be revised for certain experiments.
cfg.physics.enableElectronicNoise = 0                           # Will be revised for certain experiments.

cfg.recon.fov = 300.0
cfg.recon.sliceThickness = cfg.scanner.detectorRowSize*cfg.scanner.sid/cfg.scanner.sdd
cfg.recon.sliceCount = 1                                        # Will be revised for certain experiments.
cfg.recon.kernelType = 'standard'
cfg.recon.mu = catsim.GetMu.GetMu('water', cfg.physics.monochromatic)[0]/10 # Will be revised for certain experiments.
cfg.recon.displayWindowMin = -250             # In HU.
cfg.recon.displayWindowMax = 250              # In HU.
cfg.recon.displayWindow = cfg.recon.displayWindowMax - cfg.recon.displayWindowMin
cfg.recon.displayLevel = (cfg.recon.displayWindowMax + cfg.recon.displayWindowMin)/2
# For publication.
cfg.recon.saveImagePictureFiles = True
# For test development - comment out if desired for published figures.
# cfg.recon.displayImagePictureAxes = True
# cfg.recon.displayImagePictureTitles = True

# Top-level cfg parameters related to control of this experiment.
cfg.waitForKeypress = False             # Wait for keypress after plot display?
cfg.do_Sim = True                       # The simulation is usually run except when only varying recon parameters.
cfg.displaySimProjectionPlots = False   # Flag to display plots to screen.
cfg.do_Recon = True                     # The recon is usually run except when only varying display parameters.

##--------- Define experiment names

# The following variable sets up all relevant parameters for experiments designed to evaluate XCIST using the CatSim logo phantom.
# Uncomment the ones you want to run.

experimentNames = [
    "01_01_Ideal",
    "01_02_ViewAliasing",
    "01_03_BeamHardening",
    "01_04_Scatter",
    "01_05_Noise",
    "01_06_AllArtifacts",
]

# The current config is the base config, and will be used as the basis each time through the loop.
base = copy.deepcopy(cfg) 

for experimentIndex in range(0, len(experimentNames)):

    # Start with the "base" config.
    cfg = copy.deepcopy(base)

    # Initialize the experiment directory.
    cfg.experimentName = experimentNames[experimentIndex]
    cfg = initializeExperimentDirectory(cfg)

    # Define the specific parameters for this experiment.
    cfg = setExperimentParameters(cfg)

    # Get the title for the recon images.
    cfg.reconImageTitle = getReconImageTitle(cfg)

    # # Print the important parameters/
    # print("***")
    # print("*** Constant configuration parameters:")
    # print("* scanner.detectorPrefilter: {}".format(cfg.scanner.detectorPrefilter))
    # print("* scanner.detectorDepth: {}".format(cfg.scanner.detectorDepth))
    # print("* scanner.detectionGain: {}".format(cfg.scanner.detectionGain))
    # print("* scanner.eNoise: {}".format(cfg.scanner.eNoise))
    # print("* protocol.spectrumScaling: {}".format(cfg.protocol.spectrumScaling))
    # print("* protocol.mA: {}".format(cfg.protocol.mA))
    # print("* physics.energyCount: {}".format(cfg.physics.energyCount))
    # print("* physics.colSampleCount: {}".format(cfg.physics.colSampleCount))
    # print("* physics.rowSampleCount: {}".format(cfg.physics.rowSampleCount))
    # print("* physics.srcXSampleCount: {}".format(cfg.physics.srcXSampleCount))
    # print("* physics.srcYSampleCount: {}".format(cfg.physics.srcYSampleCount))
    # print("* physics.viewSampleCount: {}".format(cfg.physics.viewSampleCount))
    # print("* physics.energyCount: {}".format(cfg.physics.energyCount))
    # print("* recon.fov: {}".format(cfg.recon.fov))
    # print("* recon.sliceCount: {}".format(cfg.recon.sliceCount))
    # print("* recon.sliceThickness: {}".format(cfg.recon.sliceThickness))
    # print("* recon.kernelType: {}".format(cfg.recon.kernelType))
    # print("* recon.displayWindowMin: {}".format(cfg.recon.displayWindowMin))
    # print("* recon.displayWindowMax: {}".format(cfg.recon.displayWindowMax))
    # print("* recon.displayWindow: {}".format(cfg.recon.displayWindow))
    # print("* recon.displayLevel: {}".format(cfg.recon.displayLevel))
    # print("***")
    # print("*** Experiment-dependent configuration parameters:")
    # print("* scanner.detectorRowsPerMod: {}".format(cfg.scanner.detectorRowsPerMod))
    # print("* scanner.detectorRowCount: {}".format(cfg.scanner.detectorRowCount))
    # print("* scanner.detectorColOffset: {}".format(cfg.scanner.detectorColOffset))
    # print("* protocol.viewsPerRotation: {}".format(cfg.protocol.viewsPerRotation))
    # print("* protocol.viewCount: {}".format(cfg.protocol.viewCount))
    # print("* protocol.stopViewId: {}".format(cfg.protocol.stopViewId))
    # print("* protocol.maxPrep: {}".format(cfg.protocol.maxPrep))
    # print("* physics.monochromatic: {}".format(cfg.physics.monochromatic))
    # print("* physics.enableQuantumNoise: {}".format(cfg.physics.enableQuantumNoise))
    # print("* physics.enableElectronicNoise: {}".format(cfg.physics.enableElectronicNoise))
    # print("* recon.mu: {}".format(cfg.recon.mu))

    # Do the experiment.
    doExperiment(cfg)
    
