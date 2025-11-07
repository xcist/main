# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import os
import numpy as np

# -----------------------------------------------------------------------
# Aim
#   Write a CatSim view to a raw data file, in single-precision floats.
#
# Inputs
#              FileName : Name of file that contains the view data.
#             ViewIndex : The index of the current view (first view = 1).
#                         Appended to file name when sinogram is individual files for each view.
#              ViewData : The view data to be written to file.
#  SplitSinogramPerView : A flag that the sinogram file(s) represent individual views (1)
#                                                                 or one big file containing all views (0).
# Outputs
#   NONE
# -----------------------------------------------------------------------

def WriteView(FileName, ViewIndex, ViewData, SplitSinogramPerView=0):
    # Check if directory exists
    filepath = os.path.dirname(FileName)
    if filepath and not os.path.exists(filepath):
        os.makedirs(filepath)

    if SplitSinogramPerView:
        FileName = f"{FileName}.{ViewIndex}"
        mode = 'wb'
    else:
        mode = 'wb' if ViewIndex == 1 else 'ab'

    try:
        with open(FileName, mode) as FileID:
            ViewData.astype('float32').tofile(FileID)
    except IOError:
        raise IOError(f"Error: Cannot open output data file {FileName}.")