# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

# Author: Paul FitzGerald
# Date: May 10, 2022
#
# Purpose: This is an XCIST "experiment file" that is used to evaluate several aspects of XCIST simulation and recon
# using the CatSim logo phantom. The default config files are used for everything except the phantom - for that, you
# will need Phantom_CatSimLogo.cfg, which should be included a "cfg" folder where you found this file, and the files
# in the "CatSimLogo_1024" folder in XCIST's "phantoms-voxelized" repository.
#
# The CatSim logo phantom used in this experiment is a voxelized phantom. Therefore, the simulation portion of this
# evaluation only evaluates the voxelized projector. However, the recon aspects of this experiment are independent of
# the phantom/projector type.
#
# This "experiment file" performs several experiments, and each experiment requires several simulations and
# reconstructions. Results are written to the folder defined by the system environment variable XCIST_UserPath.
# Sub-folders are for each sim/recon. DON'T make this path point to a location within your local copy of the CatSim
# repository, or the output files (which are numerous and large) might end up in the repository!
#
# The experiments performed include:
# 01. Noise simulation evaluation (12 simulations/reconstructions)
# 02. Reconstruction kernal evaluation (1 simulation and 5 reconstructions from the resulting sinogram) 
# 03. Rotation speed simulation evaluation (3 simulations/reconstructions)
# 04. Number of views simulation evaluation (3 simulations/reconstructions)
# 05. Source, detector, and view sampling simulation evaluation (9 simulations/reconstructions)
# 06. X-ray scatter simulation evaluation (3 simulations/reconstructions)
# 07. Number of reconstructed slices evaluation (3 simulations/reconstructions)
# 08. Phantom offset simulation evaluation (5 simulations/reconstructions)
# 09. Reconstruction offset evaluation (1 simulation and 3 reconstructions)
#
# Each sim/recon is independent except experiments 02 and 09, which use the same sim for multiple recons.
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
import shutil
import catsim as xc
from my_commonTools import *

def setExperimentParameters(cfg):

    experimentName = cfg.experimentName

    # Some experiments use low mA.
    if experimentName == "01_01_Baseline"                        \
    or experimentName == "01_02_Physics_eNoiseOn"                \
    or experimentName == "01_03_Physics_qNoiseOn"                \
    or experimentName == "01_04_Physics_NoiseOn"                 \
    or experimentName == "01_05_Physics_1ebin"                   \
    or experimentName == "01_06_Physics_eNoiseOn_1ebin"          \
    or experimentName == "01_07_Physics_qNoiseOn_1ebin"          \
    or experimentName == "01_08_Physics_NoiseOn_1ebin"           \
    or experimentName == "01_09_Physics_Monoenergetic"           \
    or experimentName == "01_10_Physics_eNoiseOn_Monoenergetic"  \
    or experimentName == "01_11_Physics_qNoiseOn_Monoenergetic"  \
    or experimentName == "01_12_Physics_NoiseOn_Monoenergetic":
        cfg.protocol.mA = 100
    
    # Some experiments use a 128mm FOV.
    if experimentName == "02_01_Physics_NoiseOn_Recon_128mmFOV_R-LKernel"                    \
    or experimentName == "02_02_Physics_NoiseOn_Recon_128mmFOV_S-LKernel"                    \
    or experimentName == "02_03_Physics_NoiseOn_Recon_128mmFOV_SoftKernel"                   \
    or experimentName == "02_04_Physics_NoiseOn_Recon_128mmFOV_StandardKernel"               \
    or experimentName == "02_05_Physics_NoiseOn_Recon_128mmFOV_BoneKernel"                   \
    or experimentName == "05_01_Physics_SourceSampling1_Recon_128mmFOV"                      \
    or experimentName == "05_02_Physics_SourceSampling2_Recon_128mmFOV"                      \
    or experimentName == "05_03_Physics_SourceSampling3_Recon_128mmFOV"                      \
    or experimentName == "05_04_Physics_DetectorSampling1_Recon_128mmFOV"                    \
    or experimentName == "05_05_Physics_DetectorSampling2_Recon_128mmFOV"                    \
    or experimentName == "05_06_Physics_DetectorSampling3_Recon_128mmFOV":
        cfg.recon.fov = 128.0

    # Some experiments use only electronic noise.
    if experimentName == "01_02_Physics_eNoiseOn"                                             \
    or experimentName == "01_06_Physics_eNoiseOn_1ebin"                                       \
    or experimentName == "01_10_Physics_eNoiseOn_Monoenergetic":
        cfg.physics.enableElectronicNoise = 1

    # Some experiments use only quantum noise.
    if experimentName == "01_03_Physics_qNoiseOn"                                             \
    or experimentName == "01_07_Physics_qNoiseOn_1ebin"                                       \
    or experimentName == "01_11_Physics_qNoiseOn_Monoenergetic":
        cfg.physics.enableQuantumNoise = 1

    # Some experiments use electonic and quantum noise and produce unique projection data.
    if experimentName == "01_04_Physics_NoiseOn"                                             \
    or experimentName == "01_08_Physics_NoiseOn_1ebin"                                       \
    or experimentName == "01_12_Physics_NoiseOn_Monoenergetic"                               \
    or experimentName == "03_01_Physics_NoiseOn_Protocol_0p5rotation"                        \
    or experimentName == "03_02_Physics_NoiseOn_Protocol_1p0rotation"                        \
    or experimentName == "03_03_Physics_NoiseOn_Protocol_2p0rotation"                        \
    or experimentName == "04_01_Physics_NoiseOn_Protocol_100views"                           \
    or experimentName == "04_02_Physics_NoiseOn_Protocol_360views"                           \
    or experimentName == "04_03_Physics_NoiseOn_Protocol_1000views":
        cfg.physics.enableQuantumNoise = 1
        cfg.physics.enableElectronicNoise = 1

    # Some experiments use electonic and quantum noise and the projection data should be common to several recons.
    if experimentName == "02_01_Physics_NoiseOn_Recon_128mmFOV_R-LKernel"                    \
    or experimentName == "02_02_Physics_NoiseOn_Recon_128mmFOV_S-LKernel"                    \
    or experimentName == "02_03_Physics_NoiseOn_Recon_128mmFOV_SoftKernel"                   \
    or experimentName == "02_04_Physics_NoiseOn_Recon_128mmFOV_StandardKernel"               \
    or experimentName == "02_05_Physics_NoiseOn_Recon_128mmFOV_BoneKernel":
        cfg.physics.enableQuantumNoise = 1
        cfg.physics.enableElectronicNoise = 1

        # Only do the sim for the first one. Otherwise, copy the relevant projection data.
        if experimentName == "02_02_Physics_NoiseOn_Recon_128mmFOV_S-LKernel"                    \
        or experimentName == "02_03_Physics_NoiseOn_Recon_128mmFOV_SoftKernel"                   \
        or experimentName == "02_04_Physics_NoiseOn_Recon_128mmFOV_StandardKernel"               \
        or experimentName == "02_05_Physics_NoiseOn_Recon_128mmFOV_BoneKernel":
            cfg.do_Sim = False
            copyFromExperimentName = "02_01_Physics_NoiseOn_Recon_128mmFOV_R-LKernel"
            copyFromPathname = os.path.join(cfg.experimentDirectory, copyFromExperimentName,  copyFromExperimentName + ".prep")
            copyToPathname = os.path.join(cfg.experimentDirectory,  experimentName,  experimentName + ".prep")
            shutil.copy2(copyFromPathname, copyToPathname)

    # Most experiments use the baseline polyenergetic spectrum and number of energy bins,
    # but a few use different spectra.
    # We need to adjust the recon mu values for each spectrum to make water = 0 HU.
    # The adjusted mu values were determined experimentally by first reconstructing with cgf.recon.mu = 0.02,
    # measuring water HU in a cental ROI (40 pixels wide, 80 pixels high),
    # and calculating the mu required to make water = 0 HU using this formula:
    # cfg.recon.mu = 0.02 + HU(water, measured)*0.02/1000

    # Baseline:
    cfg.recon.mu = 0.019672

    # Some experiments use a polyenergetic spectrum but only one energy bin.
    if experimentName == "01_05_Physics_1ebin"                                      \
    or experimentName == "01_06_Physics_eNoiseOn_1ebin"                             \
    or experimentName == "01_07_Physics_qNoiseOn_1ebin"                             \
    or experimentName == "01_08_Physics_NoiseOn_1ebin":
        cfg.physics.energyCount = 1
        cfg.recon.mu = 0.02061

    # Some experiments use a monoenergetic spectrum.
    if experimentName == "01_01_Baseline"                                           \
    or experimentName == "01_09_Physics_Monoenergetic"                              \
    or experimentName == "01_10_Physics_eNoiseOn_Monoenergetic"                     \
    or experimentName == "01_11_Physics_qNoiseOn_Monoenergetic"                     \
    or experimentName == "01_12_Physics_NoiseOn_Monoenergetic":
        cfg.physics.monochromatic = 70
        cfg.recon.mu = 0.019326

    # Vary the recon kernel
    if experimentName == "02_01_Physics_NoiseOn_Recon_128mmFOV_R-LKernel":
        cfg.recon.kernelType = "R-L"
    if experimentName == "02_02_Physics_NoiseOn_Recon_128mmFOV_S-LKernel":
        cfg.recon.kernelType = "S-L"
    if experimentName == "02_03_Physics_NoiseOn_Recon_128mmFOV_SoftKernel":
        cfg.recon.kernelType = "Soft"
    if experimentName == "02_04_Physics_NoiseOn_Recon_128mmFOV_StandardKernel":
        cfg.recon.kernelType = "Standard"
    if experimentName == "02_05_Physics_NoiseOn_Recon_128mmFOV_BoneKernel":
        cfg.recon.kernelType = "Bone"

    # Vary the rotation time.
    if experimentName == "03_01_Physics_NoiseOn_Protocol_0p5rotation":
        cfg.protocol.rotationTime = 0.5
    if experimentName == "03_02_Physics_NoiseOn_Protocol_1p0rotation":
        cfg.protocol.rotationTime = 1.0
    if experimentName == "03_03_Physics_NoiseOn_Protocol_2p0rotation":
        cfg.protocol.rotationTime = 2.0

    # Vary the number of views.
    if experimentName == "04_01_Physics_NoiseOn_Protocol_100views":
        cfg.protocol.viewsPerRotation = 100
    if experimentName == "04_02_Physics_NoiseOn_Protocol_360views":
        cfg.protocol.viewsPerRotation = 360
    if experimentName == "04_03_Physics_NoiseOn_Protocol_1000views":
        cfg.protocol.viewsPerRotation = 1000
    cfg.protocol.viewCount = cfg.protocol.viewsPerRotation
    cfg.protocol.stopViewId = cfg.protocol.viewCount - 1

    # Vary the in-plane sampling.
    # Source sampling: sourece needs to be large to see the effect.
    if experimentName == "05_01_Physics_SourceSampling1_Recon_128mmFOV" \
    or experimentName == "05_02_Physics_SourceSampling2_Recon_128mmFOV" \
    or experimentName == "05_03_Physics_SourceSampling3_Recon_128mmFOV":
        cfg.scanner.focalspotWidth = 5.0
        cfg.scanner.focalspotLength = 5.0
    if experimentName == "05_01_Physics_SourceSampling1_Recon_128mmFOV":
        cfg.physics.srcXSampleCount = 1
    if experimentName == "05_02_Physics_SourceSampling2_Recon_128mmFOV":
        cfg.physics.srcXSampleCount = 2
    if experimentName == "05_03_Physics_SourceSampling3_Recon_128mmFOV":
        cfg.physics.srcXSampleCount = 3

    # Detector sampling: expect no effect with voxelized phantoms.
    if experimentName == "05_04_Physics_DetectorSampling1_Recon_128mmFOV":
        cfg.physics.colSampleCount = 1
    if experimentName == "05_05_Physics_DetectorSampling2_Recon_128mmFOV":
        cfg.physics.colSampleCount = 2
    if experimentName == "05_06_Physics_DetectorSampling3_Recon_128mmFOV":
        cfg.physics.colSampleCount = 3

    # View sampling: need to zoom in on a radial edge at the edge of the FOV to see the effect.
    if experimentName == "05_07_Physics_ViewSampling1_Recon_300mmFOV":
        cfg.physics.viewSampleCount = 1
    if experimentName == "05_08_Physics_ViewSampling2_Recon_300mmFOV":
        cfg.physics.viewSampleCount = 2
    if experimentName == "05_09_Physics_ViewSampling3_Recon_300mmFOV":
        cfg.physics.viewSampleCount = 3

    if experimentName == "05_01_Physics_SourceSampling1_Recon_128mmFOV"    \
    or experimentName == "05_02_Physics_SourceSampling2_Recon_128mmFOV"    \
    or experimentName == "05_03_Physics_SourceSampling3_Recon_128mmFOV"    \
    or experimentName == "05_04_Physics_DetectorSampling1_Recon_128mmFOV"  \
    or experimentName == "05_05_Physics_DetectorSampling2_Recon_128mmFOV"  \
    or experimentName == "05_06_Physics_DetectorSampling3_Recon_128mmFOV"  \
    or experimentName == "05_07_Physics_ViewSampling1_Recon_300mmFOV"      \
    or experimentName == "05_08_Physics_ViewSampling2_Recon_300mmFOV"      \
    or experimentName == "05_09_Physics_ViewSampling3_Recon_300mmFOV":
        cfg.recon.displayWindowMin = -100             # In HU.
        cfg.recon.displayWindowMax = 1300             # In HU.

    # Vary scatter
    if experimentName == "06_00_Scanner_64rows_Physics_NoScatter"      \
    or experimentName == "06_01_Scanner_64rows_Physics_ScatterScale1"  \
    or experimentName == "06_02_Scanner_64rows_Physics_ScatterScale8"  \
    or experimentName == "06_03_Scanner_64rows_Physics_ScatterScale64":
        cfg.phantom.filename = 'CatSimLogo_1024_128mmZ.json' # The phantom scale factor is 0.5, resulting in 64-mm Z.
        cfg.scanner.detectorRowsPerMod = 64
        cfg.scanner.detectorRowCount = cfg.scanner.detectorRowsPerMod
        cfg.recon.sliceCount = 64
        cfg.recon.displayWindowMin = -200             # In HU.
        cfg.recon.displayWindowMax = 200              # In HU.

    if experimentName == "06_01_Scanner_64rows_Physics_ScatterScale1"  \
    or experimentName == "06_02_Scanner_64rows_Physics_ScatterScale8"  \
    or experimentName == "06_03_Scanner_64rows_Physics_ScatterScale64":
        cfg.physics.scatterCallback = "Scatter_ConvolutionModel"
        cfg.physics.scatterKernelCallback = ""
    if experimentName == "06_01_Scanner_64rows_Physics_ScatterScale1":
        cfg.physics.scatterScaleFactor = 1
    if experimentName == "06_02_Scanner_64rows_Physics_ScatterScale8":
        cfg.physics.scatterScaleFactor = 8
    if experimentName == "06_03_Scanner_64rows_Physics_ScatterScale64":
        cfg.physics.scatterScaleFactor = 64
    
    # Some experiments use a 16-slice sim, and most of those use a 16-slice recon.
    if experimentName == "07_01_Scanner_16rows_Recon_1slice"                     \
    or experimentName == "07_02_Scanner_16rows_Recon_2slices"                    \
    or experimentName == "07_03_Scanner_16rows_Recon_16slices"                   \
    or experimentName == "08_01_Scanner_16rows_Phantom_offset0"                  \
    or experimentName == "08_02_Scanner_16rows_Phantom_offset+50mmX"             \
    or experimentName == "08_03_Scanner_16rows_Phantom_offset+50mmY"             \
    or experimentName == "08_04_Scanner_16rows_Phantom_offset+4mmZ"              \
    or experimentName == "08_05_Scanner_16rows_Phantom_offset+8mmZ"              \
    or experimentName == "09_01_Scanner_16rows_Recon_0p5mmSlices_offset0"        \
    or experimentName == "09_02_Scanner_16rows_Recon_0p5mmSlices_offset+22mmX"   \
    or experimentName == "09_03_Scanner_16rows_Recon_0p5mmSlices_offset+22mmY"   \
    or experimentName == "09_04_Scanner_16rows_Recon_0p5mmSlices_offset+1mmZ":
        cfg.scanner.detectorRowsPerMod = 16
        cfg.scanner.detectorRowCount = cfg.scanner.detectorRowsPerMod
        cfg.recon.sliceCount = 16
    # But some of them vary the number of recon slices. These only require one simulation.
    if experimentName == "07_01_Scanner_16rows_Recon_1slice":
        cfg.recon.sliceCount = 1
    if experimentName == "07_02_Scanner_16rows_Recon_2slices":
        cfg.recon.sliceCount = 2
    if experimentName == "07_03_Scanner_16rows_Recon_16slices":
        cfg.recon.sliceCount = 16

    if experimentName == "07_02_Scanner_16rows_Recon_2slices" \
    or experimentName == "07_03_Scanner_16rows_Recon_16slices":
        cfg.do_Sim = False
        copyFromExperimentName = "07_01_Scanner_16rows_Recon_1slice"
        copyFromPathname = os.path.join(cfg.experimentDirectory,  copyFromExperimentName,  copyFromExperimentName + ".prep")
        copyToPathname = os.path.join(cfg.experimentDirectory,  experimentName,  experimentName + ".prep")
        shutil.copy2(copyFromPathname, copyToPathname)

    # For phantom and recon offset tests, use 360 views, a monenergetic spectrum, and don't use oversampling.
    if experimentName == "08_01_Scanner_16rows_Phantom_offset0"               \
    or experimentName == "08_02_Scanner_16rows_Phantom_offset+50mmX"          \
    or experimentName == "08_03_Scanner_16rows_Phantom_offset+50mmY"          \
    or experimentName == "08_04_Scanner_16rows_Phantom_offset+4mmZ"           \
    or experimentName == "08_05_Scanner_16rows_Phantom_offset+8mmZ"           \
    or experimentName == "09_01_Scanner_16rows_Recon_0p5mmSlices_offset0"        \
    or experimentName == "09_02_Scanner_16rows_Recon_0p5mmSlices_offset+22mmX"   \
    or experimentName == "09_03_Scanner_16rows_Recon_0p5mmSlices_offset+22mmY"   \
    or experimentName == "09_04_Scanner_16rows_Recon_0p5mmSlices_offset+1mmZ":
        cfg.protocol.viewsPerRotation = 360
        cfg.protocol.viewCount = cfg.protocol.viewsPerRotation
        cfg.protocol.stopViewId = cfg.protocol.viewCount - 1
        cfg.physics.monochromatic = 70
        cfg.physics.srcXSampleCount = 1
        cfg.physics.srcYSampleCount = 1
        cfg.physics.rowSampleCount = 1
        cfg.physics.colSampleCount = 1
        cfg.physics.viewSampleCount = 1
        
    # Using a 16-slice sim and recon, vary the phantom offset.
    if experimentName == "08_01_Scanner_16rows_Phantom_offset0":
        cfg.phantom.centerOffset = [0.0, 0.0, 0.0]    
    if experimentName == "08_02_Scanner_16rows_Phantom_offset+50mmX":
        cfg.phantom.centerOffset = [50.0, 0.0, 0.0]    
    if experimentName == "08_03_Scanner_16rows_Phantom_offset+50mmY":
        cfg.phantom.centerOffset = [0.0, 50.0, 0.0]    
    if experimentName == "08_04_Scanner_16rows_Phantom_offset+4mmZ":
        cfg.phantom.centerOffset = [0.0, 0.0, 4.0]    
    if experimentName == "08_05_Scanner_16rows_Phantom_offset+8mmZ":
        cfg.phantom.centerOffset = [0.0, 0.0, 8.0]

    # Using a 16-slice sim and recon with 0.5-mm slices, vary the recon X, Y, and Z offset.
    if experimentName == "09_01_Scanner_16rows_Recon_0p5mmSlices_offset0"      \
    or experimentName == "09_02_Scanner_16rows_Recon_0p5mmSlices_offset+22mmX" \
    or experimentName == "09_03_Scanner_16rows_Recon_0p5mmSlices_offset+22mmY" \
    or experimentName == "09_04_Scanner_16rows_Recon_0p5mmSlices_offset+1mmZ":
        cfg.recon.sliceThickness = 0.5

    if experimentName == "09_01_Scanner_16rows_Recon_0p5mmSlices_offset0":
        cfg.recon.centerOffset = [0.0, 0.0, 0.0]    
    if experimentName == "09_02_Scanner_16rows_Recon_0p5mmSlices_offset+22mmX":
        cfg.recon.centerOffset = [22.0, 0.0, 0.0]    
    if experimentName == "09_03_Scanner_16rows_Recon_0p5mmSlices_offset+22mmY":
        cfg.recon.centerOffset = [0.0, 22.0, 0.0]
    if experimentName == "09_04_Scanner_16rows_Recon_0p5mmSlices_offset+1mmZ":
        cfg.recon.centerOffset = [0.0, 0.0, 1.0]

    # Only do the sim for the first one. Otherwise, copy the relevant projection data.
    if experimentName == "09_02_Scanner_16rows_Recon_0p5mmSlices_offset+22mmX" \
    or experimentName == "09_03_Scanner_16rows_Recon_0p5mmSlices_offset+22mmY" \
    or experimentName == "09_04_Scanner_16rows_Recon_0p5mmSlices_offset+1mmZ":
        cfg.do_Sim = False
        copyFromExperimentName = "09_01_Scanner_16rows_Recon_0p5mmSlices_offset0"
        copyFromPathname = os.path.join(cfg.experimentDirectory,  copyFromExperimentName,  copyFromExperimentName + ".prep")
        copyToPathname = os.path.join(cfg.experimentDirectory,  experimentName,  experimentName + ".prep")
        shutil.copy2(copyFromPathname, copyToPathname)

    cfg.recon.displayWindow = cfg.recon.displayWindowMax - cfg.recon.displayWindowMin
    cfg.recon.displayLevel = (cfg.recon.displayWindowMax + cfg.recon.displayWindowMin)/2

    # For single-row simulations, use the native slice thickness.
    if cfg.scanner.detectorRowCount == 1:
        cfg.recon.sliceThickness = cfg.scanner.detectorRowSize*cfg.scanner.sid/cfg.scanner.sdd
    
    return cfg


def getReconImageTitle(cfg):

    experimentName = cfg.experimentName

    # Most experiments have only one image, so don't add that to the title.
    cfg.addSliceInfoToReconImageTitle = False

    if experimentName == "01_01_Baseline"                                          \
    or experimentName == "01_02_Physics_eNoiseOn"                                  \
    or experimentName == "01_03_Physics_qNoiseOn"                                  \
    or experimentName == "01_04_Physics_NoiseOn"                                   \
    or experimentName == "01_05_Physics_1ebin"                                     \
    or experimentName == "01_06_Physics_eNoiseOn_1ebin"                            \
    or experimentName == "01_07_Physics_qNoiseOn_1ebin"                            \
    or experimentName == "01_08_Physics_NoiseOn_1ebin"                             \
    or experimentName == "01_09_Physics_Monoenergetic"                             \
    or experimentName == "01_10_Physics_eNoiseOn_Monoenergetic"                    \
    or experimentName == "01_11_Physics_qNoiseOn_Monoenergetic"                    \
    or experimentName == "01_12_Physics_NoiseOn_Monoenergetic":
        formatString = WLString(cfg) + "cfg.physics.monochromatic = {}; cfg.physics.energyCount = {};"
        string1 = formatString.format(cfg.physics.monochromatic, cfg.physics.energyCount)
        formatString = "cfg.physics.enableElectronicNoise = {}; cfg.protocol.spectrumScaling = {}"
        string2 = formatString.format(cfg.physics.enableElectronicNoise, cfg.protocol.spectrumScaling)
        formatString = "cfg.physics.enableQuantumNoise = {}; cfg.protocol.mA = {}; cfg.recon.mu = {}"
        string3 = formatString.format(cfg.physics.enableQuantumNoise, cfg.protocol.mA, cfg.recon.mu)
        return string1 + "\n" + string2 + "\n" + string3

    if experimentName == "02_01_Physics_NoiseOn_Recon_128mmFOV_R-LKernel"           \
    or experimentName == "02_02_Physics_NoiseOn_Recon_128mmFOV_S-LKernel"           \
    or experimentName == "02_03_Physics_NoiseOn_Recon_128mmFOV_SoftKernel"          \
    or experimentName == "02_04_Physics_NoiseOn_Recon_128mmFOV_StandardKernel"      \
    or experimentName == "02_05_Physics_NoiseOn_Recon_128mmFOV_BoneKernel":
        formatString = WLString(cfg) + "cfg.recon.kernelType = {:s}"
        return formatString.format(cfg.recon.kernelType)

    if experimentName == "03_01_Physics_NoiseOn_Protocol_0p5rotation"               \
    or experimentName == "03_02_Physics_NoiseOn_Protocol_1p0rotation"               \
    or experimentName == "03_03_Physics_NoiseOn_Protocol_2p0rotation":
        formatString = WLString(cfg) + "cfg.protocol.rotationTime = {} s"
        return formatString.format(cfg.protocol.rotationTime)

    if experimentName == "04_01_Physics_NoiseOn_Protocol_100views"                  \
    or experimentName == "04_02_Physics_NoiseOn_Protocol_360views"                  \
    or experimentName == "04_03_Physics_NoiseOn_Protocol_1000views":
        formatString = WLString(cfg) + "cfg.protocol.viewsPerRotation = {}"
        return formatString.format(cfg.protocol.viewsPerRotation)

    if experimentName == "05_01_Physics_SourceSampling1_Recon_128mmFOV"             \
    or experimentName == "05_02_Physics_SourceSampling2_Recon_128mmFOV"             \
    or experimentName == "05_03_Physics_SourceSampling3_Recon_128mmFOV":
        formatString = WLString(cfg) +  "cfg.physics.srcXSampleCount = {}"
        return formatString.format(cfg.physics.srcXSampleCount)

    if experimentName == "05_04_Physics_DetectorSampling1_Recon_128mmFOV"           \
    or experimentName == "05_05_Physics_DetectorSampling2_Recon_128mmFOV"           \
    or experimentName == "05_06_Physics_DetectorSampling3_Recon_128mmFOV":
        formatString = WLString(cfg) +  "cfg.physics.colSampleCount = {}"
        return formatString.format(cfg.physics.colSampleCount)

    if experimentName == "05_07_Physics_ViewSampling1_Recon_300mmFOV"               \
    or experimentName == "05_08_Physics_ViewSampling2_Recon_300mmFOV"               \
    or experimentName == "05_09_Physics_ViewSampling3_Recon_300mmFOV":
        formatString = WLString(cfg) +  "cfg.physics.viewSampleCount = {}"
        return formatString.format(cfg.physics.viewSampleCount)

    if experimentName == "06_00_Scanner_64rows_Physics_NoScatter"           \
    or experimentName == "06_01_Scanner_64rows_Physics_ScatterScale1"       \
    or experimentName == "06_02_Scanner_64rows_Physics_ScatterScale8"       \
    or experimentName == "06_03_Scanner_64rows_Physics_ScatterScale64":
        # For multi-slice scans, add slice info to the title.
        cfg.addSliceInfoToReconImageTitle = True
        formatString = WLString(cfg) + "cfg.physics.scatterScaleFactor = {}"
        return formatString.format(cfg.physics.scatterScaleFactor)

    if experimentName == "07_01_Scanner_16rows_Recon_1slice"      \
    or experimentName == "07_02_Scanner_16rows_Recon_2slices"     \
    or experimentName == "07_03_Scanner_16rows_Recon_16slices":
        # For multi-slice scans, add slice info to the title.
        cfg.addSliceInfoToReconImageTitle = True
        formatString = "cfg.scanner.detectorRowCount = {}; cfg.recon.sliceCount = {}"
        return formatString.format(cfg.scanner.detectorRowCount, cfg.recon.sliceCount)

    if experimentName == "08_01_Scanner_16rows_Phantom_offset0"      \
    or experimentName == "08_02_Scanner_16rows_Phantom_offset+50mmX" \
    or experimentName == "08_03_Scanner_16rows_Phantom_offset+50mmY" \
    or experimentName == "08_04_Scanner_16rows_Phantom_offset+4mmZ"  \
    or experimentName == "08_05_Scanner_16rows_Phantom_offset+8mmZ":
        # For multi-slice scans, add slice info to the title.
        cfg.addSliceInfoToReconImageTitle = True
        formatString = "cfg.phantom.centerOffset[X, Y, Z] = [{}, {}, {}]"
        return formatString.format(cfg.phantom.centerOffset[0], cfg.phantom.centerOffset[1], cfg.phantom.centerOffset[2])

    if experimentName == "09_01_Scanner_16rows_Recon_0p5mmSlices_offset0"      \
    or experimentName == "09_02_Scanner_16rows_Recon_0p5mmSlices_offset+22mmX" \
    or experimentName == "09_03_Scanner_16rows_Recon_0p5mmSlices_offset+22mmY" \
    or experimentName == "09_04_Scanner_16rows_Recon_0p5mmSlices_offset+1mmZ":
        # For multi-slice scans, add slice info to the title.
        cfg.addSliceInfoToReconImageTitle = True
        formatString = "cfg.recon.centerOffset[X, Y, Z] = [{}, {}, {}]"
        return formatString.format(cfg.recon.centerOffset[0], cfg.recon.centerOffset[1], cfg.recon.centerOffset[2])


##--------- Initialize

userPath = getUserPath()
xc.CommonTools.my_path.add_search_path(userPath)

# Use the default cfg parameters found in the default .cfg files, except use a specific phantom file.

cfg = xc.CatSim(xc.CommonTools.my_path.find("cfg", "Phantom_CatSimLogo.cfg", ""))
cfg.experimentDirectory = os.path.join(userPath, "examples", "evaluation", "experiment_01_CatSimLogo_SimAndRecon")

# These are changes to the defaults config parameters to be used for the "base" experiment,
# but some of these and some other defaults are overridden for certain experiments.

cfg.scanner.detectorRowsPerMod = 1                              # Will be revised for certain experiments.
cfg.scanner.detectorRowCount = cfg.scanner.detectorRowsPerMod   # Will be revised for certain experiments.
cfg.scanner.detectorColOffset = 0.25

cfg.protocol.viewsPerRotation = 1000                            # Will be revised for certain experiments.
cfg.protocol.viewCount = cfg.protocol.viewsPerRotation          # Will be revised for certain experiments.
cfg.protocol.stopViewId = cfg.protocol.viewCount - 1            # Will be revised for certain experiments.
cfg.protocol.spectrumScaling = 0.931
cfg.protocol.mA = 300                                           # Will be revised for certain experiments.

cfg.physics.enableQuantumNoise = 0                              # Will be revised for certain experiments.
cfg.physics.enableElectronicNoise = 0                           # Will be revised for certain experiments.
cfg.physics.energyCount = 12 # Using 120 kVp, so 10 kV/bin.

cfg.recon.fov = 300.0                                           # Will be revised for certain experiments.
cfg.recon.sliceCount = 1                                        # Will be revised for certain experiments.
cfg.recon.kernelType = 'Standard'                               # Will be revised for certain experiments.
cfg.recon.mu = 0.02                                             # Will be revised for certain experiments.
cfg.recon.displayWindowMin = -200             # In HU.
cfg.recon.displayWindowMax = 200              # In HU.
cfg.recon.displayWindow = cfg.recon.displayWindowMax - cfg.recon.displayWindowMin
cfg.recon.displayLevel = (cfg.recon.displayWindowMax + cfg.recon.displayWindowMin)/2
# For publication.
cfg.recon.saveImagePictureFiles = True
# For test development - comment out if desired for published figures.
# cfg.recon.displayImagePictureAxes = True
# cfg.recon.displayImagePictureTitles = True
# For single-row simulations, use the native slice thickness.     Will be revised for certain experiments.
if cfg.scanner.detectorRowCount == 1:
    cfg.recon.sliceThickness = cfg.scanner.detectorRowSize*cfg.scanner.sid/cfg.scanner.sdd

# Top-level cfg parameters related to control of this experiment.
cfg.waitForKeypress = False             # Wait for keypress after plot display?
cfg.do_Sim = True                       # The simulation is usually run except when only varying recon parameters.
cfg.displaySimProjectionPlots = False   # Flag to display plots to screen.
cfg.do_Recon = True                     # The recon is usually run except when only varying display parameters.

##--------- Define experiment names

# The following variable sets up all relevant parameters for experiments designed to evaluate XCIST using the CatSim logo phantom.
# Uncomment the ones you want to run.

experimentNames = [
    # "01_01_Baseline",
    # "01_02_Physics_eNoiseOn",
    # "01_03_Physics_qNoiseOn",
    # "01_04_Physics_NoiseOn",
    # "01_05_Physics_1ebin",
    # "01_06_Physics_eNoiseOn_1ebin",
    # "01_07_Physics_qNoiseOn_1ebin",
    # "01_08_Physics_NoiseOn_1ebin",
    # "01_09_Physics_Monoenergetic",
    # "01_10_Physics_eNoiseOn_Monoenergetic",
    # "01_11_Physics_qNoiseOn_Monoenergetic",
    # "01_12_Physics_NoiseOn_Monoenergetic",

    # "02_01_Physics_NoiseOn_Recon_128mmFOV_R-LKernel", # Needs to be done before the next 4 because those use projections from this.
    # "02_02_Physics_NoiseOn_Recon_128mmFOV_S-LKernel",
    # "02_03_Physics_NoiseOn_Recon_128mmFOV_SoftKernel",
    # "02_04_Physics_NoiseOn_Recon_128mmFOV_StandardKernel",
    # "02_05_Physics_NoiseOn_Recon_128mmFOV_BoneKernel",

    # "03_01_Physics_NoiseOn_Protocol_0p5rotation",
    # "03_02_Physics_NoiseOn_Protocol_1p0rotation",
    # "03_03_Physics_NoiseOn_Protocol_2p0rotation",
    
    # "04_01_Physics_NoiseOn_Protocol_100views",
    # "04_02_Physics_NoiseOn_Protocol_360views",
    # "04_03_Physics_NoiseOn_Protocol_1000views",
    
    # "05_01_Physics_SourceSampling1_Recon_128mmFOV",
    # "05_02_Physics_SourceSampling2_Recon_128mmFOV",
    # "05_03_Physics_SourceSampling3_Recon_128mmFOV",
    # "05_04_Physics_DetectorSampling1_Recon_128mmFOV",
    # "05_05_Physics_DetectorSampling2_Recon_128mmFOV",
    # "05_06_Physics_DetectorSampling3_Recon_128mmFOV",
    # "05_07_Physics_ViewSampling1_Recon_300mmFOV",
    # "05_08_Physics_ViewSampling2_Recon_300mmFOV",
    # "05_09_Physics_ViewSampling3_Recon_300mmFOV",

    "06_00_Scanner_64rows_Physics_NoScatter",
    "06_01_Scanner_64rows_Physics_ScatterScale1",
    "06_02_Scanner_64rows_Physics_ScatterScale8",
    "06_03_Scanner_64rows_Physics_ScatterScale64",

    "07_01_Scanner_16rows_Recon_1slice",
    "07_02_Scanner_16rows_Recon_2slices",
    "07_03_Scanner_16rows_Recon_16slices",

    "08_01_Scanner_16rows_Phantom_offset0",
    "08_02_Scanner_16rows_Phantom_offset+50mmX",
    "08_03_Scanner_16rows_Phantom_offset+50mmY",
    "08_04_Scanner_16rows_Phantom_offset+4mmZ",
    "08_05_Scanner_16rows_Phantom_offset+8mmZ",

    "09_01_Scanner_16rows_Recon_0p5mmSlices_offset0",  # Needs to be done before the next 3 because those use projections from this.
    "09_02_Scanner_16rows_Recon_0p5mmSlices_offset+22mmX",
    "09_03_Scanner_16rows_Recon_0p5mmSlices_offset+22mmY",
    "09_04_Scanner_16rows_Recon_0p5mmSlices_offset+1mmZ",
]

# The current config is the base config, and will be used as the basis each time through the loop.
base = copy.deepcopy(cfg) 

for experimentIndex in range(0, len(experimentNames)):

    # Start with the "base" config.
    cfg = copy.deepcopy(base)

    cfg.experimentName = experimentNames[experimentIndex]

    # Initialize the experiment directory.
    cfg = initializeExperimentDirectory(cfg)

    # Define the specific parameters for this experiment.
    cfg = setExperimentParameters(cfg)

    # Get the title for the recon images.
    cfg.reconImageTitle = getReconImageTitle(cfg)

    # Do the experiment.
    doExperiment(cfg)
    
