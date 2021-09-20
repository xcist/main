######################################################################################################################
# dicom_to_voxelized_phantom
# This file reads in the folder pathname for a series of DICOM images and a list of materials and generates a 
# voxelized phantom based on those images and materials. The mu value of each material is calculated via geant
# and the DICOM image is separated with each voxel identified with one particular material. Smaller differences
# are treated as different densities of the same material (i.e. dense bone or osteoporosis)
# Depends on installation of catsim, and currently the directory for CatSim is hard coded.
# This must be fixed before distribution
# Also - the current file does not receive the necessary data. That needs to be corrected as well
# Author: Chad Bircher 8/18/2021
#
#
# Needed Updates: 
# Build CatSim - should include the update to the source, so try removing that line
# load settings from CFG file
# Include a definition of how the phantom was defined from the CFG file
# Allow the CFG file to pass threshold values for material decomposition
######################################################################################################################
from pathlib import Path
from os import listdir
from os.path import isfile, join, getctime
import numpy as np
import re
import pydicom
import copy
import json
import matplotlib.pyplot as plt # possibly update to a prettier ploting tool such as plot.ly - TBD
######################################################################################################################
# This was placed before I set up CatSim correctly on my machine. It should not generally be needed
# import sys
# sys.path.append('C:\\Users\\212312408\\Documents\\GitHub\\CatSim') 
######################################################################################################################
from catsim.GetMu  import GetMu
from catsim.CommonTools import source_cfg

class IndexTracker:
    # Tracker for allowing scroll wheel to move through a large number of images
    # The images are assumed to be stored as [image_number, x, y]
    # In the event that the images are stored as [x, y, image number] update [self.ind, :, :] -> [:, :, self.ind]
    def __init__(self, ax, X):
        self.ax = ax
        self.X = X
        self.slices, rows, cols = X.shape
        self.ind = self.slices//2
        self.im = ax.imshow(self.X[self.ind, :, :]) # initialize the image number
        self.update() # draw the default image (typically the central slice)

    def on_scroll(self, event): # when the wheel has scrolled
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices # increment image number by one
        else:
            self.ind = (self.ind - 1) % self.slices # decrement image number by one
        self.update() # update the image

    def update(self): # When updating
        self.im.set_data(self.X[self.ind, :, :]) # select the image
        self.ax.set_ylabel('slice %s' % self.ind) # title the y axis with the image number
        self.im.axes.figure.canvas.draw() # update the figure with the selected image

def generate_vp(dicom_example, materials, num_slices, phtName):
    num_materials = len(materials)
    x_pixels = len(dicom_example.pixel_array) # number of x pixels
    y_pixels = len(dicom_example.pixel_array[0]) # number of y pixels
    wrapper_file = {
        "construction_description": 'Some description here, possibly the cfg file??',
        "n_materials": num_materials,  # phantom made of 2 materials, so there're 2 density volume files
        "mat_name": materials['names'],  # material names
        "mu_values": materials['mu_values'],  # Geant defined mu values 
        "mu_thresholds": materials['threshold_values'], # Thresholds 
        "volumefractionmap_filename": materials['file_names'],  # filenames of density volumes
        "volumefractionmap_datatype":["float"] * num_materials,  # binary data type
        "cols":[x_pixels] * num_materials,  # dimensions
        "rows":[y_pixels] * num_materials,
        "slices":[num_slices] * num_materials,  # Total number of slices: one slice for each image for each material
        "x_size":[dicom_example.PixelSpacing[0]] * num_materials,  # voxel size: X direction
        "y_size":[dicom_example.PixelSpacing[1]] * num_materials,  # voxel size: Y direction
        "z_size":[dicom_example.SliceThickness] * num_materials, #If there are failures here try SpacingBetweenSlices or extract from SeriesDescription
        "x_offset":[0.5 + dicom_example.Columns / 2] * num_materials,  # center position, in voxels
        "y_offset":[0.5 + dicom_example.Rows / 2] * num_materials,
        "z_offset":[0.5 + num_slices / 2] * num_materials
    }
    with open(phtName + '.json', 'w') as outfile:
        json.dump(wrapper_file, outfile, indent=4)  # print out the json wrapper file with pretty printing
    return wrapper_file


def define_materials(phantom, tmp, num_slices):
    # Materials definition section
    ee = 63.948  # water Mu = 0.2 at this energy - hard coded in generate_phantom
    mu_list0 = []
    density_file = {}
    fnames = []
    materialList = phantom.materials
    phtName = phantom.phantomname
    for index, material in enumerate(materialList):
        fname = phtName + '.density_' + str(index)  # name for this particular material in the phantom group
        fnames.append(fname)                    # This can probably be combined with the line above - add the name to the list of names
        mu_list0.append(GetMu(material, ee)[0]) # append the list of mu with the next material
        density_file.update({material: np.zeros((num_slices, len(tmp.pixel_array), len(tmp.pixel_array[0])), dtype=np.float32)})  # add a placeholder for the phantom of this material
        # Caluclate the Thresholds - mu boundaries between materials
        if material == materialList[0]:
            thld_list = [0]
        else:
            thld_list.append(mu_list0[index-1]*0.55 + mu_list0[index]*0.45) # 55% lower mu material, 45% upper mu material
    # sort the mu_list0 and keep the index of the new values
    indices = sorted(range(len(mu_list0)),key=mu_list0.__getitem__)
    # for each index calculate the threshold between the previous and current materials
    # match the order for the mu, materials, and thresholds
    mu_list1 = [mu_list0[index] for index in indices]
    materialList1 = [materialList[index] for index in indices]
    for index, material in enumerate(materialList1):
        if material == materialList1[0]:
            thld_list1 = [0]
        else:
            thld_list1.append(mu_list1[index-1]*0.55 + mu_list1[index]*0.45) # 55% lower mu material, 45% upper mu material

    material_dict = {'names': materialList, 'file_names': fnames, 'mu_values': mu_list1, 'threshold_values': thld_list1} # material dictionary: perhaps add the mu list and threshold list to this dictionary?
    return density_file, material_dict


def phantom_run_decision(mypath, phantom):
    # read in files in current folder
    # check if the phantom for density_0 exists, if so check it's creation time
    run_phantom = True
    phantom_exists = False
    files0 = listdir()
    if phantom.phantomname + '.density_0' in files0:
        time_phantom = getctime(phantom.phantomname + '.density_0')
        run_phantom = False
        phantom_exists = True

    # Find all files and DICOM files using the assigned directory
    allfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    dcmfiles = [j for j in allfiles if '.dcm' in j]

    # Sort DICOM files: Initial order is alphabetical - switch to numeric
    # Also check whether to run the DICOM phantom: 
    #   run if there is no DICOM phantom or if the phantom is older than the DICOM images
    #   do not run if the phantom is newer than the DICOM images
    #   Print an error message if there is no phantom and there are no DICOM images

    dcm_num = []
    tmp = []
    if len(dcmfiles) > 0:
        # read in one dcm file to use as a sample for definitions
        fname = mypath + '\\' + dcmfiles[0]
        tmp = pydicom.dcmread(fname)
        for fname in dcmfiles:
            dcm_num.append(int(re.findall(r'\d+', fname)[0])) # dicom number in the file name - alphabetical order - convert to number
            time_dicom = getctime(mypath + '/' + fname)
        indices = sorted(range(len(dcm_num)),key=dcm_num.__getitem__)  # sort the dicom numbers by numeric order
        new_dcm = [dcmfiles[index] for index in indices]  # find the index associated with the new sorted number to open files in numeric rather than alphabetic order
        try:
            # check if ranges for slices are given
            try:
                 test = list(range(phantom.sliceRange[0],1 + phantom.sliceRange[1])) # If a single range is given use that range of slices
            except:
                test = []
                for j in phantom.sliceRange:
                    test += list(range(j[0], 1+j[1])) # if a series of ranges are given then include the slices in each range
            new_dcm = [new_dcm[i] for i in test] # filter the dicom data on the selected range(s)
        except:
            new_dcm = new_dcm # keep the original if the ranges are inappropriate or if sliceRanges are not given

        if phantom_exists:
            if time_phantom < time_dicom: # Phantom exists, but the DICOM data have been updated - update phantom
                # print('New DICOM data, calculating Phantom')
                run_phantom = True
            else:
                # print('Using existing phantom files')
                run_phantom = False
        else: # DICOM files exist but there is no matching phantom. Build phantom
            # print('No phantom exists, building one using DICOM data')
            run_phantom = True
    else:
        if phantom_exists:
            # print('No DICOM data available, using existing phantom file')
            pass
        else:
            # print('No phantom file exists and no DICOM data are available')
            pass
    
    # Given current code design considerations always run the phantom. Retaining code above in case our priorities change.
    run_phantom = True

    return new_dcm, run_phantom, tmp

def generate_phantom(new_dcm, mypath, material_dict, density_file):
    thld_list = material_dict['threshold_values']
    materialList = material_dict['names']
    mu_list0 = material_dict['mu_values']
    # Generate Phantom Arrays: 
    #   Read in each DICOM file
    #   Separate into materials using thresholds calculated above
    #   Calculate relative density
    #   Append to arrays linked to materials
    mu_water = 0.2
    for dcm_index, partial_name in enumerate(new_dcm): # has size limitations on share drive (532MB) - always run locally
        fname = mypath + '\\' + partial_name
        tmp = pydicom.dcmread(fname)                    # Read DICOM File
        tmp_array = tmp.pixel_array * mu_water / 1000   # Calculate density relative to water
        tmp_array[tmp_array < 0] = 0                    # Remove air
        bounds = copy.deepcopy(thld_list)               # Needed to prevent the next line from appending to thld_list
        bounds.append(1.1*tmp_array.max())              # Set the upper bound for the last material to include the highest pixel in the array
        for material_index, material in enumerate(materialList):
            array0 = copy.deepcopy(tmp_array)                          # Read in the full array for this slice
            array0[array0 < bounds[material_index]] = 0         # Zero out pixels below lower threshold
            array0[array0 > bounds[material_index+1]] = 0       # Zero out pixels above upper threshold
            array0 *= 1/mu_list0[material_index]                # Calculate density of material in phantom
            density_file[material][dcm_index] = array0     # Store the density in the density dictionary
    return density_file

def write_files(material_dict, new_dcm, phtName, tmp, density_file):
    fnames = material_dict['file_names']
    materialList = material_dict['names']
    # Write Phantom Files
    for index, fname in enumerate(fnames):
        with open(fname, 'wb') as fout:
            fout.write(density_file[materialList[index]])  # Write the phantom layer to the matching file name
    # Generate and print wrapper file
    generate_vp(tmp, material_dict, len(new_dcm), phtName)
    return 0

def ConvertDICOM(phantom):
    phtName = phantom.phantomname  # Name of the DICOM data folder, also the name of the phantom attached to that folder
    mypath = phantom.dicomdirectory + phtName # parent folder for the DICOM data folder - is this located correctly?
    # Check for existing phantom and DICOM files
    new_dcm, run_phantom, sample_dicom = phantom_run_decision(mypath, phantom)
    # new_dcm is the list of DICOM files in numeric order
    # run_phantom is a boolean used to decide whether to calculate the new phantom
    # sample_dicom is one of the DICOM files: the information in that file is used to:
    #   determine the size of the placeholder arrays in define_material
    #   define many of the variables printed to the wrapper file in write_files

    # If the phantom does not exist or if the phantom is older than the DICOM images (re)run the phantom
    if run_phantom:
        density_file, material_dict = define_materials(phantom, sample_dicom, len(new_dcm))
        # check if thresholds have been entered into the config file. If so check that there are the correct number of thresholds. If so use them.
        try: 
            if len(phantom.thresholds) == len(material_dict['threshold_values']): # If user has defined thresholds and the number is correct
                material_dict.update({'threshold_values': phantom.thresholds})    # Then update the material dictionary with user defined thresholds
        except:
            print('using calculated thresholds')                                  # Otherwise use the calculated thresholds
        density_file = generate_phantom(new_dcm, mypath, material_dict, density_file)    # Generate the density maps for the phantoms
        write_files(material_dict, new_dcm, phtName, sample_dicom, density_file)         # Write the data (both phantoms and wrapper) to the appropriate files

        try:
            if phantom.showPhantom:                           # If the user has defined showPhantom as True
                num_materials = len(material_dict['names'])   # find the number of materials
                if  num_materials <= 3:                       # With 3 or fewer materials define a 1 by x grid
                    rows = 1
                    cols = num_materials
                elif num_materials == 4:                      # With 4 materials define a 2x2 grid
                    rows = 2
                    cols = 2
                elif num_materials <= 6:                      # With 5 or 6 materials define a 2x3 grid
                    rows = 2
                    cols = 3
                # Here I should add grids for 7-8 (2x4), 9-12 (3x4), 13-16 (4x4), and 17-25(5x5) materials. Beyond 25 is very unlikely to be useful or meaningful
                fig, ax = plt.subplots(rows, cols, sharex=True, sharey=True) # Define a plot with rows and axes defined above, all axes are linked for zooming purposes
                tracker = []
                for plot_num, material in enumerate(material_dict['names']):  # For each material
                    # Identify the subplot
                    if num_materials <= 3:
                        this_axis = ax[plot_num%cols]
                    else:
                        col = plot_num%cols
                        this_axis = ax[(num_materials - col) / rows, col]

                    tracker.append(IndexTracker(this_axis,density_file[material])) # Link this subplot to the tracker
                    fig.canvas.mpl_connect('scroll_event', tracker[plot_num].on_scroll)  # define event of scrolling to change the view slice
                    this_axis.set_title(material) # name the subplot as the material name
                plt.show()
        except:
            pass
    return 0

if __name__ == "__main__":
    # read in config file
    test_config = source_cfg('C:\\Users\\212312408\\Documents\\GitHub\\CatSim\\catsim\\cfg\\Phantom_Head.cfg')
    # When integrating into CatSim the location of the phantom file needs to be passed with the call statement
    ConvertDICOM(test_config.phantom)
