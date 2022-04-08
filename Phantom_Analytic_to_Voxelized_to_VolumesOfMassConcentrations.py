# Generated with SMOP  0.41
from libsmop import *
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m

    # -----------------------------------------------------------------------
#   Program Name: Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m
#   Authors:  Samit Basu, Bruno De Man, and Jed Pack (GE Global Research)
#   Organization: GE Global Research and GE Healthcare, General Electric Company
#   Version:  6.0.3
#   Date:  Feb 3, 2015
#   Class: GE Confidential. General Electric Proprietary Data (c) 2012 General Electric Company
    
    # Aim
#   This script runs catvoxel if needed, to produce a set of material volumes, and invokes dose_map_v. 
# 
# Input
#   Since this is a script (not a function) all variables visible in the code that invokes this are also visible here.
#   Assumed: A valid cfg stucture, includng cfg.phantom_filename containing a valid pathname to an analytic phantom.
    
    # Output
#   A voxelized phantom file, equivilent to the analytic phantom described in cfg.phantom_filename at entry.
#   cfg.phantom_filename is updated with the name of the voxelized phantom.
#  
# History:
#   2012-09-20 Paul FitzGerald (GE Global Research)
#              Renamed - was dose_map_p.
#              Cleaned up and added "Verbose" output. 
#              Changed force_translation to cfg.force_phantom_conversion.
#   2013-04-17 Paul FitzGerald (GE Global Research)
#              Fixed cfg.force_phantom_conversion variable name bug.
# -----------------------------------------------------------------------
    
    global Verbose
    Directory,FileName,Extension=fileparts(cfg.phantom_filename,nargout=3)
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:31
    if logical_not((strcmp(Extension,'.pp') or strcmp(Extension,'.ppm'))):
        error('Unknown phantom file format: %s',cfg.phantom_filename)
    
    VoxelizedPhantomFilename=strrep(cfg.phantom_filename,Extension,'.vp')
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:36
    if logical_not(exist(VoxelizedPhantomFilename)) or first_file_newer(cfg.phantom_filename,VoxelizedPhantomFilename) or (isfield(cfg,'force_phantom_conversion') and cfg.force_phantom_conversion):
        Headline(cellarray(['',' ','Voxelizing analytic phantom...',' ']),Verbose.Phantom,- 1)
        Headline(cellarray([' ','Setting phantom matrix per cfg.phantom_samples_*',' ']),Verbose.Phantom)
        adjust='cfg.material_volumes = 1;'
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:45
        adjust=concat([adjust,' cfg.Nx   = cfg.phantom_samples_xy;'])
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:46
        adjust=concat([adjust,' cfg.dx   = cfg.phantom_samples_voxelsize;'])
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:47
        adjust=concat([adjust,' cfg.xoff = (cfg.Nx+1)/2;'])
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:48
        adjust=concat([adjust,' cfg.Ny   = cfg.phantom_samples_xy;'])
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:49
        adjust=concat([adjust,' cfg.dy   = cfg.phantom_samples_voxelsize;'])
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:50
        adjust=concat([adjust,' cfg.yoff = (cfg.Ny+1)/2;'])
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:51
        adjust=concat([adjust,' cfg.Nz   = cfg.phantom_samples_z;'])
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:52
        adjust=concat([adjust,' cfg.dz   = cfg.phantom_samples_voxelsize;'])
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:53
        adjust=concat([adjust,' cfg.zoff = (cfg.Nz+1)/2;'])
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:54
        cfg.write_vp = copy(1)
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:56
        Volume,MaterialList=catvoxel([],cfg,adjust,logical_not(Verbose.Phantom),nargout=2)
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:57
        Headline(cellarray([' ','... done voxelizing analytic phantom.',' ','']),Verbose.Phantom,1)
    
    cfg.phantom_filename = copy(VoxelizedPhantomFilename)
# Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.m:63
    Phantom_Voxelized_to_VolumesOfMassConcentrations