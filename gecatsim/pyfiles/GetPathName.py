# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

# -----------------------------------------------------------------------
# Aim:
#   This function returns the full path name to a specified file name, if that file is on the Python path.
#   If not, an empty string is returned.
# -----------------------------------------------------------------------

import os
import sys

def GetPathName(FileName):
    Paths = sys.path
    for Path in Paths:
        PathName = os.path.join(Path, FileName)
        if os.path.exists(PathName):
            return PathName
    return ''