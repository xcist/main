# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

# Author: Paul FitzGerald
# Date: April 18, 2022
#
# Purpose: These are common functions used b several experiment files.

import os
import numpy as np
import gecatsim as xc
import gecatsim.reconstruction.pyfiles.recon as recon

def initializeExperimentDirectory(cfg):

    resultsPath = os.path.join(cfg.experimentDirectory, cfg.experimentName)
    
    print("**************************")
    print("* Experiment: {:s}".format(cfg.experimentName))
    print("**************************")

    # Create the results folder if it doesn't exist.
    if not os.path.exists(resultsPath):
        os.makedirs(resultsPath)
        
    cfg.resultsName = os.path.join(resultsPath, cfg.experimentName)

    return cfg


def doExperiment(cfg):

    # Run the simulation.
    runSim(cfg)

    # Run the recon.
    cfg = runRecon(cfg)

    return cfg


def WLString(cfg):
    # Used below to create strings with the window/Level info.

    if cfg.recon.unit == 'HU':
        return "W/L = {:.0f}/{:.0f} {:s}; ".format(cfg.recon.displayWindow, cfg.recon.displayLevel, cfg.recon.unit)
    if cfg.recon.unit == '/cm':
        return "W/L = {:.2f}/{:.2f} {:s}; ".format(cfg.recon.displayWindow, cfg.recon.displayLevel, cfg.recon.unit)
    if cfg.recon.unit == '/mm':
        return "W/L = {:.3f}/{:.3f} {:s}; ".format(cfg.recon.displayWindow, cfg.recon.displayLevel, cfg.recon.unit)


def runSim(cfg):

    import matplotlib.pyplot as plt

    if cfg.do_Sim == True:
        print("**************************")
        print("* Running the simulation *")
        print("**************************")

        cfg.run_all()  # as specified by protocol.scanTypes

    ##--------- Show selected results - 4 views for each of selected rows.

    if cfg.displaySimProjectionPlots:
    
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
            # Plot the first 3 rows, a middle row, and the last 3 rows.
            rowIndicesToPlot = [0, 1, 2, int(rowCount/2), rowCount-3, rowCount-2, rowCount-1]

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

        experimentName = cfg.experimentName
        if rowCount == 16:
            if experimentName == "08_01_Scanner_16rows_Phantom_offset0"       \
            or experimentName == "08_02_Scanner_16rows_Phantom_offset+50mmX"  \
            or experimentName == "08_03_Scanner_16rows_Phantom_offset+50mmY":
                # Confirm that symmetrical rows are exactly the same.
                rowsToDiff = np.array([[0, rowCount-1], [1, rowCount-2], [2, rowCount-3]])
            elif experimentName == "08_04_Scanner_16rows_Phantom_offset+4mmZ":
                # Compare the 4 middle rows.
                firstRowToDiff = int(rowCount/2 - 1)
                rowsToDiff = np.array([[firstRowToDiff,   firstRowToDiff+1],
                                    [firstRowToDiff+1, firstRowToDiff+2],
                                    [firstRowToDiff+2, firstRowToDiff+3]])
            elif experimentName == "08_05_Scanner_16rows_Phantom_offset+8mmZ":
                # Compare the last 5 rows.
                rowsToDiff = np.array([[rowCount-5, rowCount-4],
                                    [rowCount-4, rowCount-3],
                                    [rowCount-3, rowCount-2],
                                    [rowCount-2, rowCount-1]])

            if experimentName == "08_01_Scanner_16rows_Phantom_offset0"       \
            or experimentName == "08_02_Scanner_16rows_Phantom_offset+50mmX"  \
            or experimentName == "08_03_Scanner_16rows_Phantom_offset+50mmY"  \
            or experimentName == "08_04_Scanner_16rows_Phantom_offset+4mmZ"   \
            or experimentName == "08_05_Scanner_16rows_Phantom_offset+8mmZ":
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

    if userPath is None:
        raise Exception("******** Error! Please set the environment variable 'XCIST_UserPath'")
    # Convert to a list to see if more than one path was specified.
    userPath = userPath.split(";")
    if len(userPath) > 1:
        raise Exception("******** Error! Environment variable 'XCIST_UserPath' can only contain a single path.\nIt contains {:s}.".format(userPath))

    # Convert back to a simple string.
    userPath = userPath[0]
    if not os.path.exists(userPath):
        raise Exception("******** Error! Environment variable 'XCIST_UserPath' not found.".format(userPath))

    return userPath
