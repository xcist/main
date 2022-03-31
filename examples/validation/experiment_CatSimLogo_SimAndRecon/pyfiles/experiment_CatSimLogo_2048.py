# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import os
import copy
import shutil
import catsim as xc
import reconstruction.pyfiles.recon as recon

def runSim(cfg):

    import numpy as np
    import matplotlib.pyplot as plt

    if cfg.do_Sim == True:
        print("**************************")
        print("* Running the simulation *")
        print("**************************")
        print("cfg.physics.enableElectronicNoise = {}".format(cfg.physics.enableElectronicNoise))
        print("cfg.physics.enableQuantumNoise = {}".format(cfg.physics.enableQuantumNoise))

        cfg.run_all()  # as specified by protocol.scanTypes

    ##--------- Show selected results - 4 views for each of selected rows.

    print("*********************************")
    print("* Plotting selected projections *")
    print("*********************************")

    projectionData = xc.rawread(cfg.resultsName + ".prep",
                    [cfg.protocol.viewCount, cfg.scanner.detectorRowCount, cfg.scanner.detectorColCount],
                    'float')

    rowCount = cfg.scanner.detectorRowCount
    if rowCount == 1:
        rowIndicesToPlot = [0]
    elif rowCount <= 8:
        # rowIndicesToPlot = range(8, 16)
        rowIndicesToPlot = range(0, rowCount)
    else:
        if (experimentName == "Phantom_offset_0"
        or experimentName == "Phantom_offset_+50-mm_X"
        or experimentName == "Phantom_offset_+50-mm_Y"):
            # Plot the first 3 rows, a middle row, and the last 3 rows.
            rowIndicesToPlot = [0, 1, 2, int(rowCount/2), rowCount-3, rowCount-2, rowCount-1]
        elif experimentName == "Phantom_offset_+4-mm_Z":
            # Plot the 2 middle rows.
            rowIndicesToPlot = [int(rowCount/2)-1, int(rowCount/2)]
        elif experimentName == "Phantom_offset_+8-mm_Z":
            # Plot the last 4 rows.
            rowIndicesToPlot = [rowCount-4, rowCount-3, rowCount-2, rowCount-1]
        else:
            rowIndicesToPlot = range(0, rowCount)

    for rowIndexToPlot in rowIndicesToPlot:
        for viewIndexToPlot in range(0, cfg.protocol.viewCount, int(cfg.protocol.viewCount/4)):
            viewToPlot = projectionData[viewIndexToPlot, rowIndexToPlot, :]
            # viewToPlot = np.subtract(projectionData[viewIndexToPlot, rowIndexToPlot, :], projectionData[viewIndexToPlot, 15, :])
            plt.figure(int(viewIndexToPlot+1))
            plt.plot(viewToPlot)
            plt.title("All " + str(cfg.scanner.detectorColCount) + " columns"
                    + " of row(s) " + str([s + 1 for s in rowIndicesToPlot]) + " of " + str(cfg.scanner.detectorRowCount)
                    + ", view " + str(viewIndexToPlot+1) + " of " + str(cfg.protocol.viewCount))
            fileName = cfg.resultsName + "_" + "View_" + str(viewIndexToPlot+1)
            plt.savefig(fileName, bbox_inches='tight')

    if rowCount == 16:
        if (experimentName == "Phantom_offset_0"
            or experimentName == "Phantom_offset_+50-mm_X"
            or experimentName == "Phantom_offset_+50-mm_Y"):
            # Confirm that symmetrical rows are exactly the same.
            rowsToDiff = np.array([[0, rowCount-1], [1, rowCount-2], [2, rowCount-3]])
        elif experimentName == "Phantom_offset_+4-mm_Z":
            # Compare the 4 middle rows.
            firstRowToDiff = int(rowCount/2 - 1)
            rowsToDiff = np.array([[firstRowToDiff,   firstRowToDiff+1],
                                [firstRowToDiff+1, firstRowToDiff+2],
                                [firstRowToDiff+2, firstRowToDiff+3]])
        elif experimentName == "Phantom_offset_+8-mm_Z":
            # Compare the last 5 rows.
            rowsToDiff = np.array([[rowCount-5, rowCount-4],
                                [rowCount-4, rowCount-3],
                                [rowCount-3, rowCount-2],
                                [rowCount-2, rowCount-1]])

            numRowsToDiff = np.shape(rowsToDiff)[0]
            for rowIndex in range(0, numRowsToDiff):
                diff = np.subtract(projectionData[0, rowsToDiff[rowIndex, 0], :], projectionData[0, rowsToDiff[rowIndex, 1], :])
                print("Max diff of rows ", rowsToDiff[rowIndex, 0], " and ", rowsToDiff[rowIndex, 1], " is ", np.max(np.abs(diff)))

    plt.draw()
    plt.pause(1)
    if cfg.waitForKeypress:
        print("*******************************************")
        print("* Press Enter to close plots and continue *")
        input("*******************************************")
    
    plt.close('all')


def runRecon(cfg):

    print("******************************")
    print("* Running the reconstruction *")
    print("******************************")

    cfg = recon.recon(cfg)

    return cfg


def getUserPath():

    # Get the user-specified environment variable.
    userPath = os.environ.get('XCIST_UserPath')

    # Convert to a list to see if more than one path was specified.
    userPath = userPath.split(";")
    if len(userPath) > 1:
        raise Exception("******** Error! Environment variable 'XCIST_UserPath' can only contain a single path.\nIt contains {:s}.".format(userPath))

    # Convert back to a simple string.
    userPath = userPath[0]
    if not os.path.exists(userPath):
        raise Exception("******** Error! Environment variable 'XCIST_UserPath' not found.".format(userPath))

    return userPath


##--------- Initialize

userPath = getUserPath()
experimentDirectory = os.path.join(userPath, "my_experiments", "experiment_CatSimLogo_SimAndRecon")

cfg = xc.CatSim(os.path.join(experimentDirectory, "cfg", "Phantom_CatSimLogo_2048.cfg"),
                os.path.join(experimentDirectory, "cfg", "Recon_CatSimLogo.cfg"))

# These are changes to the defaults config parameters to be used for the "base" experiment,
# but some of these and some other defaults are overridden for certain experiments.

cfg.scanner.detectorRowsPerMod = 1
cfg.scanner.detectorRowCount = cfg.scanner.detectorRowsPerMod
cfg.scanner.detectorColOffset = 0.25

cfg.physics.enableQuantumNoise = 0
cfg.physics.enableElectronicNoise = 0
cfg.physics.energyCount = 12 # Using 120 kVp so 10 kV/bin.

cfg.protocol.viewsPerRotation = 1000 # About the minimum to get minimal aliasing
cfg.protocol.viewCount = cfg.protocol.viewsPerRotation
cfg.protocol.stopViewId = cfg.protocol.viewCount-1
cfg.protocol.spectrumScaling = 0.931
cfg.protocol.mA = 300

cfg.recon.fov = 300.0
cfg.recon.displayImagePictures = True
cfg.recon.saveImagePictureFiles = True
cfg.recon.saveImageVolume = True
cfg.recon.saveSingleImages = True

# The following variable sets up all relevant parameters for experiments designed to evaluate XCIST geometry using the CatSim logo phantom.
# Uncomment the ones you want to run.

experimentNames = [
    "CatSimLogo_2048",
    # "Physics_eNoiseOn",
    # "Physics_qNoiseOn",
    # "Physics_NoiseOn",
    # "Physics_1ebin",
    # "Physics_eNoiseOn_1ebin",
    # "Physics_qNoiseOn_1ebin",
    # "Physics_NoiseOn_1ebin",
    # "Physics_Monoenergetic",
    # "Physics_eNoiseOn_Monoenergetic",
    # "Physics_qNoiseOn_Monoenergetic",
    # "Physics_NoiseOn_Monoenergetic",

    # "Physics_NoiseOn_Recon_128mmFOV_R-LKernel", # Needs to be done before the next 4 because those use projections from this.
    # "Physics_NoiseOn_Recon_128mmFOV_S-LKernel",
    # "Physics_NoiseOn_Recon_128mmFOV_SoftKernel",
    # "Physics_NoiseOn_Recon_128mmFOV_StandardKernel",
    # "Physics_NoiseOn_Recon_128mmFOV_BoneKernel",

    # "Physics_NoiseOn_Protocol_0p5rotation",
    # "Physics_NoiseOn_Protocol_1p0rotation",
    # "Physics_NoiseOn_Protocol_2p0rotation",
    
    # "Physics_NoiseOn_Protocol_100views",
    # "Physics_NoiseOn_Protocol_360views",
    # "Physics_NoiseOn_Protocol_1000views",
    
    # "Physics_SourceSampling1_Recon_128mmFOV",
    # "Physics_SourceSampling2_Recon_128mmFOV",
    # "Physics_SourceSampling3_Recon_128mmFOV",
    # "Physics_DetectorSampling1_Recon_128mmFOV",
    # "Physics_DetectorSampling2_Recon_128mmFOV",
    # "Physics_DetectorSampling3_Recon_128mmFOV",
    # "Physics_ViewSampling1_Recon_300mmFOV",
    # "Physics_ViewSampling2_Recon_300mmFOV",
    # "Physics_ViewSampling3_Recon_300mmFOV",

    # "Physics_ScatterScale0p5",
    # "Physics_ScatterScale1p0",
    # "Physics_ScatterScale2p0",

    # "Scanner_16slices_Recon_1slices",
    # "Scanner_16slices_Recon_2slices",
    # "Scanner_16slices_Recon_16slice",

    # "16slices_Phantom_offset+50-mm_X",
    # "16slices_Phantom_offset+50-mm_Y",
    # "16slices_Phantom_offset+4-mm_Z",
    # "16slices_Phantom_offset+8-mm_Z",
]

# The current config is the base config, and will be used as the basis each time through the loop.
base = copy.deepcopy(cfg) 


for experimentIndex in range(0, len(experimentNames)):

    cfg = copy.deepcopy(base)
    cfg.waitForKeypress = False             # Wait for keypress after plot display?
    cfg.do_Sim = False                       # The simulation is usually run except when only varying recon parameters.
    # cfg.displayWindowMin = -1000             # In HU.
    # cfg.displayWindowMax = -997              # In HU.


    print("**************************")
    print("* Base config *")
    print("**************************")
    print("cfg.physics.enableElectronicNoise = {}".format(cfg.physics.enableElectronicNoise))
    print("cfg.physics.enableQuantumNoise = {}".format(cfg.physics.enableQuantumNoise))

    experimentName = experimentNames[experimentIndex]
    resultsPath = os.path.join(experimentDirectory, experimentName)

    # Create the results folder if it doesn't exist.
    if not os.path.exists(resultsPath):
        os.makedirs(resultsPath)
    cfg.resultsName = os.path.join(resultsPath, experimentName)

    # Some experiments use a 128-mm FOV.
    if experimentName == "Physics_NoiseOn_Recon_128mmFOV_R-LKernel"                    \
    or experimentName == "Physics_NoiseOn_Recon_128mmFOV_S-LKernel"                    \
    or experimentName == "Physics_NoiseOn_Recon_128mmFOV_SoftKernel"                   \
    or experimentName == "Physics_NoiseOn_Recon_128mmFOV_StandardKernel"               \
    or experimentName == "Physics_NoiseOn_Recon_128mmFOV_BoneKernel"                   \
    or experimentName == "Physics_SourceSampling1_Recon_128mmFOV"                      \
    or experimentName == "Physics_SourceSampling3_Recon_128mmFOV"                      \
    or experimentName == "Physics_DetectorSampling1_Recon_128mmFOV"                    \
    or experimentName == "Physics_DetectorSampling3_Recon_128mmFOV":
        cfg.recon.fov = 128.0

    # Some experiments use only electronic noise.
    if experimentName == "Physics_eNoiseOn"                                             \
    or experimentName == "Physics_eNoiseOn_1ebin"                                       \
    or experimentName == "Physics_eNoiseOn_Monoenergetic":
        cfg.physics.enableElectronicNoise = 1

    # Some experiments use only quantum noise.
    if experimentName == "Physics_qNoiseOn"                                             \
    or experimentName == "Physics_qNoiseOn_1ebin"                                       \
    or experimentName == "Physics_qNoiseOn_Monoenergetic":
        cfg.physics.enableQuantumNoise = 1

    # Some experiments use electonic and quantum noise and produce unique projection data.
    if experimentName == "Physics_NoiseOn"                                             \
    or experimentName == "Physics_NoiseOn_1ebin"                                       \
    or experimentName == "Physics_NoiseOn_Monoenergetic"                               \
    or experimentName == "Physics_NoiseOn_Protocol_0p5rotation"                        \
    or experimentName == "Physics_NoiseOn_Protocol_1p0rotation"                        \
    or experimentName == "Physics_NoiseOn_Protocol_2p0rotation"                        \
    or experimentName == "Physics_NoiseOn_Protocol_100views"                           \
    or experimentName == "Physics_NoiseOn_Protocol_360views"                           \
    or experimentName == "Physics_NoiseOn_Protocol_1000views":
        cfg.physics.enableQuantumNoise = 1
        cfg.physics.enableElectronicNoise = 1

    # Some experiments use electonic and quantum noise and the projection data should be common to several recons.
    if experimentName == "Physics_NoiseOn_Recon_128mmFOV_R-LKernel"                    \
    or experimentName == "Physics_NoiseOn_Recon_128mmFOV_S-LKernel"                    \
    or experimentName == "Physics_NoiseOn_Recon_128mmFOV_SoftKernel"                   \
    or experimentName == "Physics_NoiseOn_Recon_128mmFOV_StandardKernel"               \
    or experimentName == "Physics_NoiseOn_Recon_128mmFOV_BoneKernel":
        cfg.physics.enableQuantumNoise = 1
        cfg.physics.enableElectronicNoise = 1

        # Only do the sim for the first one. Otherwise, copy the relevant projection data.
        if experimentName == "Physics_NoiseOn_Recon_128mmFOV_S-LKernel"                    \
        or experimentName == "Physics_NoiseOn_Recon_128mmFOV_SoftKernel"                   \
        or experimentName == "Physics_NoiseOn_Recon_128mmFOV_StandardKernel"               \
        or experimentName == "Physics_NoiseOn_Recon_128mmFOV_BoneKernel":
            copyFromExperimentName = "Physics_NoiseOn_Recon_128mmFOV_R-LKernel"
            cfg.do_Sim = False
    
        if not cfg.do_Sim:
            copyFromPathname = os.path.join(cfg.resultsPath,  copyFromExperimentName + ".prep")
            copyToPathname = cfg.resultsName + ".prep"
            shutil.copy2(copyFromPathname, copyToPathname)

    # Some experiments use a polyenergetic spectrum but only one energy bin.
    if experimentName == "Physics_1ebin"                                       \
    or experimentName == "Physics_eNoiseOn_1ebin"                              \
    or experimentName == "Physics_qNoiseOn_1ebin"                              \
    or experimentName == "Physics_NoiseOn_1ebin":
        cfg.physics.energyCount = 1

    # Some experiments use a monoenergetic spectrum.
    if experimentName == "Physics_Monoenergetic"                                       \
    or experimentName == "Physics_eNoiseOn_Monoenergetic"                              \
    or experimentName == "Physics_qNoiseOn_Monoenergetic"                              \
    or experimentName == "Physics_NoiseOn_Monoenergetic":
        cfg.physics.monochromatic = 70
        cfg.protocol.spectrumScaling = 1.0

    # Vary the recon kernel
    if experimentName == "Physics_NoiseOn_Recon_128mmFOV_R-LKernel":
        recon.kernelType = "R-L"
    if experimentName == "Physics_NoiseOn_Recon_128mmFOV_S-LKernel":
        cfg.recon.kernelType = "S-L"
    if experimentName == "Physics_NoiseOn_Recon_128mmFOV_SoftKernel":
        cfg.recon.kernelType = "Soft"
    if experimentName == "Physics_NoiseOn_Recon_128mmFOV_StandardKernel":
        cfg.recon.kernelType = "Standard"
    if experimentName == "Physics_NoiseOn_Recon_128mmFOV_BoneKernel":
        cfg.recon.kernelType = "Bone"

    # Vary rotation time
    if experimentName == "Physics_NoiseOn_Protocol_0p5rotation":
        cfg.protocol.rotationTime = 0.5
    if experimentName == "Physics_NoiseOn_Protocol_1p0rotation":
        cfg.protocol.rotationTime = 1.0
    if experimentName == "Physics_NoiseOn_Protocol_2p0rotation":
        cfg.protocol.rotationTime = 2.0

    # Vary number of views
    if experimentName == "Physics_NoiseOn_Protocol_100views":
        cfg.protocol.viewsPerRotation = 100
    if experimentName == "Physics_NoiseOn_Protocol_360views":
        cfg.protocol.viewsPerRotation = 360
    if experimentName == "Physics_NoiseOn_Protocol_1000views":
        cfg.protocol.viewsPerRotation = 1000
    cfg.protocol.viewCount = cfg.protocol.viewsPerRotation
    cfg.protocol.stopViewId = cfg.protocol.viewCount-1

    # Vary in-plane sampling
    if experimentName == "Physics_SourceSampling1_Recon_128mmFOV":
        cfg.physics.srcXSampleCount = 1
    if experimentName == "Physics_SourceSampling2_Recon_128mmFOV":
        cfg.physics.srcXSampleCount = 2
    if experimentName == "Physics_SourceSampling3_Recon_128mmFOV":
        cfg.physics.srcXSampleCount = 3
    if experimentName == "Physics_DetectorSampling1_Recon_128mmFOV":
        cfg.physics.rowSampleCount = 1
    if experimentName == "Physics_DetectorSampling2_Recon_128mmFOV":
        cfg.physics.rowSampleCount = 2
    if experimentName == "Physics_DetectorSampling3_Recon_128mmFOV":
        cfg.physics.rowSampleCount = 3
    if experimentName == "Physics_ViewSampling1_Recon_300mmFOV":
        cfg.physics.viewSampleCount = 1
    if experimentName == "Physics_ViewSampling2_Recon_300mmFOV":
        cfg.physics.viewSampleCount = 2
    if experimentName == "Physics_ViewSampling3_Recon_300mmFOV":
        cfg.physics.viewSampleCount = 3

    # Vary scatter
    if experimentName == "Physics_ScatterScale0p5" \
    or experimentName == "Physics_ScatterScale1p0" \
    or experimentName == "Physics_ScatterScale2p0":
        cfg.physics.scatterCallback = "Scatter_ConvolutionModel"
        cfg.physics.scatterKernelCallback = "scatter_kernel.dat"
    if experimentName == "Physics_ScatterScale0p5":
        cfg.physics.scatterScaleFactor = 0.5
    if experimentName == "Physics_ScatterScale1p0":
        cfg.physics.scatterScaleFactor = 1.0
    if experimentName == "Physics_ScatterScale2p0":
        cfg.physics.scatterScaleFactor = 2.0
    
    # Some experiments use a 16-slice sim, and most of those use a 16-slice recon.
    if experimentName == "Scanner_16slices_Recon_1slices"  \
    or experimentName == "Scanner_16slices_Recon_2slices"  \
    or experimentName == "Scanner_16slices_Recon_16slice"  \
    or experimentName == "16slices_Phantom_offset0"        \
    or experimentName == "16slices_Phantom_offset+50-mm_X" \
    or experimentName == "16slices_Phantom_offset+50-mm_Y" \
    or experimentName == "16slices_Phantom_offset+4-mm_Z"  \
    or experimentName == "16slices_Phantom_offset+8-mm_Z":
        cfg.scanner.detectorRowsPerMod = 16
        cfg.scanner.detectorRowCount = cfg.scanner.detectorRowsPerMod
        cfg.recon.sliceCount = 16
    # But some of them vary the number of recon slices. These only require one simulation.
    if experimentName == "Scanner_16slices_Recon_1slice":
        cfg.recon.sliceCount = 1
        copyFromExperimentName = experimentName
    if experimentName == "Scanner_16slices_Recon_2slices":
        cfg.recon.sliceCount = 2
    if experimentName == "Scanner_16slices_Recon_16slices":
        cfg.recon.sliceCount = 16

    if experimentName == "Scanner_16slices_Recon_2slices" \
    or experimentName == "Scanner_16slices_Recon_16slices":
        copyFromPathname = os.path.join(experimentDirectory,  copyFromExperimentName,  copyFromExperimentName + ".prep")
        copyToPathname = cfg.resultsName + ".prep"
        shutil.copy2(copyFromPathname, copyToPathname)

    # Using a 16-slice sim and recon, vary the phantom offset
    if experimentName == "Phantom_Offset0":
        cfg.phantom.centerOffset = [0.0, 0.0, 0.0]    
    if experimentName == "Phantom_Offset+50mmX":
        cfg.phantom.centerOffset = [50.0, 0.0, 0.0]    
    if experimentName == "Phantom_Offset+50mmY":
        cfg.phantom.centerOffset = [0.0, 50.0, 0.0]    
    if experimentName == "Phantom_Offset+4mmZ":
        cfg.phantom.centerOffset = [0.0, 0.0, 4.0]    
    if experimentName == "Phantom_Offset+8mmZ":
        cfg.phantom.centerOffset = [0.0, 0.0, 8.0]

    print("**************************")
    print("* Experiment: {:s}".format(experimentName))
    print("**************************")
    print("**************************")
    print("* Final config *")
    print("**************************")
    print("cfg.physics.enableElectronicNoise = {}".format(cfg.physics.enableElectronicNoise))
    print("cfg.physics.enableQuantumNoise = {}".format(cfg.physics.enableQuantumNoise))

    runSim(cfg)
    cfg = runRecon(cfg)


