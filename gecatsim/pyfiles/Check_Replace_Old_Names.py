'''
This script checks for old variable names and replaces them with new variable names.
The original source file is modified in place.

By default it will do this search-and-replace in the current directory.
This can be changed via the current_working_dir parameter below.
'''

from pathlib import Path, PosixPath
import glob
import sys

this_file = sys.argv[0]

# Name mapping - ["NEW_NAME", "OLD_NAME"]
new_old_names_list = [ 
    ("phantom.definitionCallback", "phantom.callback"),
    ("phantom.projectorThreadCount", "phantom.projectorNumThreads"),
    ("physics.detColSampleCount", "physics.colSampleCount"),
    ("physics.detRowSampleCount", "physics.rowSampleCount"),
    ("physics.recalcRayAngles", "physics.recalcRayAngle"),
    ("physics.recalcPhant", "physics.recalcPht"),
    ("physics.xrayCrosstalkCallback", "physics.crosstalkCallback"),
    ("protocol.startViewNr", "protocol.startViewId"),
    ("protocol.stopViewNr", "protocol.stopViewId"),
    ("protocol.spectrumScaleCurrent", "protocol.spectrumUnit_mm"),
    ("protocol.spectrumScaleArea", "protocol.spectrumUnit_mA"),
    ("protocol.bowtieFilename", "protocol.bowtie"),
    ("protocol.maxP", "protocol.maxPrep"),
    ("recon.fieldOfView", "recon.fov"),
    ("recon.matrixSize", "recon.imageSize"),
    ("recon.algorithm", "recon.reconType"),
    ("recon.kernel", "recon.kernelType"),
    ("recon.startAngle", "recon.unit"),
    ("recon.HUMu", "recon.mu"),
    ("recon.HUOffset", "recon.huOffset")
]

# Dictionary of new and old names
new_old_names_dict = dict(new_old_names_list)

new_names = new_old_names_dict.keys()
old_names = new_old_names_dict.values()

# Directory of files where the variables needs to be replaced. Change as needed
current_working_dir = "."

# File extensions of files where the variables need to be replaced.
file_types = ['.py', '.cfg', '.cpp']

# Get all the files and directories under the current working directory
matching_files = [matching_file for matching_file in Path(current_working_dir).rglob('*') if matching_file.suffix in file_types]

# Take care not to replace items in this file
matching_files.remove(PosixPath('check_replace_old_names.py'))

# Consider only those files with the above listed extensions
for matching_file in matching_files:

    print(f"\nCurrent file : {matching_file}")

    file_changed = False
    changed_file_content = ''

    with open(matching_file, 'r') as inFile:
        for index, current_line in enumerate(inFile):
            replaced = False
            # Check all the old variable names in the current line
            for old_name in old_names:
                if old_name in current_line:
                    file_changed = True
                    matching_new_name = list(new_old_names_dict.keys())[list(new_old_names_dict.values()).index(old_name)]
                    # Replace the old variable name with new name
                    changed_file_content = changed_file_content + current_line.replace(old_name, matching_new_name)
                    replaced = True
                    print(f"Line {index}. Changed old name '{old_name}' to  new name = {matching_new_name}")
            # If not replaced with any new variables then take the original line of file
            if not replaced:
                changed_file_content = changed_file_content + current_line

    # Save the file if any variable name is replaced
    if file_changed:
        file = open(matching_file, "w")
        file.write(changed_file_content)

