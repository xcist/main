# Generated with SMOP  0.41
# Phantom_Analytic_Get.m

    # -----------------------------------------------------------------------
#   Program Name: Phantom_Analytic_Get.m                                   
#   Authors:  Samit Basu, Bruno De Man, Jed Pack (GE Global Research)
#   Organization: GE Global Research and GE Healthcare, General Electric Company
#   Version:  6.0.3
#   Date:  Feb 3, 2015
#   Class: GE Confidential. General Electric Proprietary Data (c) 2012 General Electric Company
    
    # Aim
#   Read and parse a phantom from file, store it in C memory, and extract the materials
    
    # Inputs
#   cfg.phantom_filename     The phantom file in one of the following formats:
#      .pp format            The original FORBILD format used by CatSim 1.0 and 2.0
#      .ppm format           A FreeMat script that defines the phantom as a cell array of structures
#                            For example, object{1}.type = 1; object{1}.half_axes = [2 1 5]; ...
#                            All fields must be specified for each object:
#                              (type, center, half_axes, euler_angs, density, shape, axial_lims, material)
#                            The file must also contain a material list such as this:
#                              materialList = {'water' 'lung' 'bone'};
    
    #   cfg.phantom_scale        A scale factor to resize phantom,
#                            also used to convert from cm to mm (cfg.phantom_scale=10).
    
    #   cfg.time                 (optional) only used for certain dynamic phantoms
    
    # Outputs
#   materialList             Cell array of phantom materials (strings)
#                            Note that this is set when the phantom file is "sourced".
    
    # History: 
#   2012-10-09 Paul FitzGerald (GE Global Research)
#              Renamed - was ParsePhantomZZ.
#              Cleaned up and added "Verbose" output.
# -----------------------------------------------------------------------

import os
import sys
import numpy as np
from catsim.parse_analytical_ppm import parse_analytical_ppm
from catsim.Phantom_Analytic_SetObjects import Phantom_Analytic_SetObjects
from catsim.Phantom_Analytic_ConvertTori import Phantom_Analytic_ConvertTori
from catsim.Phantom_Analytic_BoundObjects import Phantom_Analytic_BoundObjects 


#NOTE: phobject here is object in matlab, phantobject here is phobject in matlab
def Phantom_Analytic_Get(cfg):

    #breakpoint()
    PhantomFilename=cfg.phantom.filename
    Scale=cfg.phantom.scale
    print('Reading phantom file {}...'.format(PhantomFilename))

    #NOTE by jiayong: following two should be a list of dictionaries
    #phantObject = []
    
    PhantomPath = os.path.dirname(PhantomFilename)
    BaseName = os.path.basename(PhantomFilename)
    PhantomBaseName, PhantomExtension = BaseName.split('.')
    if 'pp' == PhantomExtension:
        ppmPhantomFilename = PhantomFilename.replace('.pp','.ppm')
        #if cfg.force_phantom_conversion or logical_not(first_file_newer(ppmPhantomFilename,PhantomFilename)):
        if cfg.force_phantom_conversion or os.path.getmtime(PhantomFilename) > os.path.getmtime(ppmPhantomFilename):
            if os.path.exist(PhantomFilename):
                Phantom_Analytic_pp_to_ppm(PhantomFilename.split('.')[0])
            else:
                sys.exit('Phantom file {} not found'.format(PhantomFilename))
    else:
        if 'ppm' == (PhantomExtension):
            ppmPhantomFilename = PhantomFilename
    
    print('Reading phantom file {}.'.format(ppmPhantomFilename))
   
    #TODO: how to source file in xcist?, ask mingye how to do this part
    # phobject is just a temporary var, its values will be saved to objs and passed to setobject
    phobject = parse_analytical_ppm(ppmPhantomFilename)
  
    #NOTE: object comes from .ppm file
    #TODO: forget for now, not used in currrent cfg
    #if isfield(cfg,'phantom_rotation_y') and any(cfg.phantom_rotation_y) != 0:
    #    cfg.phantom_rotation_y_analytic = copy(cfg.phantom_rotation_y)
    #
    #if isfield(cfg,'phantom_rotation_y_analytic') and cfg.phantom_rotation_y_analytic != 0:
    #    ang=cfg.phantom_rotation_y_analytic
    #    sina=sin(dot(- ang / 180,pi))
    #    cosa=cos(dot(- ang / 180,pi))
    #    transmat=concat([[cosa,0,sina],[0,1,0],[- sina,0,cosa]])
    #    for ii in arange(1,length(object.type)).reshape(-1):
    #        if any(object.euler_angs(ii,arange(2,3)) != 0) or all(object.euler_angs(ii,1) != concat([0,90,- 90])):
    #            error('CatSim: phantom rotation currently only works without original rotation')
    #        object.center[ii,arange()]=dot(object.center(ii,arange()),transmat)
    #        object.euler_angs[ii,arange()]=concat([90,ang,90])
    #        if logical_not(isempty(object.clip[ii])):
    #            for jj in arange(1,size(object.clip[ii],1)).reshape(-1):
    #                if all(abs(object.clip[ii](jj,arange(1,3))) == concat([1,0,0])):
    #                    object.clip[ii][jj,arange(1,3)]=dot(object.clip[ii](jj,arange(1,3)),transmat)
    #
    ## phantom position offset, Mingye
    #if isfield(cfg,'phantom_position_offset') and any(cfg.phantom_position_offset) != 0:
    #    cfg.phantom_position_offset_analytic = copy(cfg.phantom_position_offset)
    #
    #if isset('cfg.phantom_position_offset_analytic') and any(cfg.phantom_position_offset_analytic) != 0:
    #    cfg.phantom_position_offset_analytic = copy(ravel(cfg.phantom_position_offset_analytic).T)
    #    for ii in arange(1,length(object.type)).reshape(-1):
    #        object.center[ii,arange()]=object.center(ii,arange()) + cfg.phantom_position_offset_analytic
    #        if logical_not(isempty(object.clip[ii])):
    #            # Mingye Wu, Sept 15 2017
    #            # get the delta D of clip distance
    #            for jj in arange(1,size(object.clip[ii],1)).reshape(-1):
    #                d=get_clip_dD(object.clip[ii](jj,arange(1,3)),cfg.phantom_position_offset_analytic)
    #                object.clip[ii][jj,4]=object.clip[ii](jj,4) + d
    if len(phobject['materialList']) < max(phobject['material'])[0]:
        sys.exit('Phantom cannot reference a material that is not in the list')
    objs = {'params': [[] for _ in range(len(phobject['type']))],
            'clip': [[] for _ in range(len(phobject['type']))]}
    if np.any(np.array(phobject['type'])[:,0]==100):
        #TODO: not do this for now
        #Phantom_Polygonal(object,Scale)
        pass
    else:
        # convert boxes to cylinders with clipping
        print('Converting boxes to cylinders with clipping.')
        for i in range(len(phobject['type'])):
            if phobject['type'][i] == 8:
                phobject['type'][i] = 2
                Bt = Rmat(phobject['euler_angs'][i]).T
                phobject['clip'][i] = np.vstack([phobject['clip'][i]],[Bt[:,0].T, phobject['half_axes'][i,0] + phobject['center'][i],Bt[:,0]])
                phobject['clip'][i] = np.vstack([phobject['clip'][i]],[-Bt[:,0].T, phobject['half_axes'][i,0] - phobject['center'][i],Bt[:,0]])
                phobject['clip'][i] = np.vstack([phobject['clip'][i]],[Bt[:,1].T, phobject['half_axes'][i,1] + phobject['center'][i],Bt[:,1]])
                phobject['clip'][i] = np.vstack([phobject['clip'][i]],[-Bt[:,1].T, phobject['half_axes'][i,1] - phobject['center'][i],Bt[:,1]])
                phobject['half_axes'][i,:2]= phobject['half_axes'][i,:2]*np.sqrt(2)

        # Scale and convert to a compact format
        print('Scaling and compacting format.')
        indsToScale= list(range(6)) + [14,15]

        for i in range(len(phobject['type'])):
            if np.logical_not(np.any(phobject['type'][i][0] in [3,7])):
                tmp = phobject['center'][i] + phobject['half_axes'][i] + phobject['euler_angs'][i] + phobject['density'][i] + phobject['type'][i] + [0] + phobject['shape'][i] + phobject['axial_lims'][i] + phobject['material'][i]
                tmp = np.array(tmp)
                tmp[indsToScale] = tmp[indsToScale]*Scale
                #breakpoint()
                objs['params'][i]=tmp
                tmp=phobject['clip'][i]
                #if np.logical_not(np.isempty(tmp)):
                if len(tmp)!=0:
                    tmp[:,3]= tmp[:,3]*Scale
                objs['clip'][i]=tmp
            else:
                objs['params'][i]= np.vstack([phobject['center'][i],phobject['half_axes'][i],phobject['euler_angs'][i],phobject['density'][i],phobject['type'][i],0,phobject['shape'][i],phobject['axial_lims'][i],phobject['material'][i]])
                objs['clip'][i]=phobject['clip'][i]
        # Set the objects, and pass to C.
        objs['params'] = np.array(objs['params'])
        #objs['clip'] = np.array(objs['clip'])
        phantObject = Phantom_Analytic_SetObjects(objs)
        objs['params'] = Phantom_Analytic_ConvertTori(phantObject,objs['params'],objs['clip'])
        phantObject = Phantom_Analytic_BoundObjects(phantObject,objs['params'])
        print('... done with phantom.')
    
    return phantObject, phobject
    
    
def Rmat(p=None):
    p=p*pi / 180
    p1=p[1]
    p2=p[2]
    p3=p[3]
    R=np.vstack([np.cos(p1)*np.cos(p3) - np.sin(p1)*np.cos(p2)*np.sin(p3), np.cos(p1)*np.sin(p3) + np.sin(p1)*np.cos(p2)*np.cos(p3), np.sin(p1)*np.sin(p2)],
            [-np.sin(p1)*np.cos(p3) - np.cos(p1)*np.cos(p2)*np.sin(p3), -np.sin(p1)*np.sin(p3) + np.cos(p1)*np.cos(p2)*np.cos(p3),np.cos(p1)*np.sin(p2)],
            [np.sin(p2)*np.sin(p3), -np.sin(p2)*np.cos(p3), np.cos(p2)])
    return R
    
# not used for now    
#def get_clip_dD(c0=None,c1=None,*args,**kwargs):
#    varargin = get_clip_dD.varargin
#    nargin = get_clip_dD.nargin
#
#    # Mingye Wu, Sept 15 2017
## get the delta D of clip
## c0: original clip vector
## c1: phantom off-centering coordinates
#    
#    r0=norm(c0,2)
#    r1=norm(c1,2)
#    d=norm(c0 - c1,2)
#    cosA=(r0 ** 2 + r1 ** 2 - d ** 2) / (dot(dot(2,r0),r1))
#    
#    cosA[logical_or(isnan(cosA),cosA) == logical_or(inf,cosA) == - inf]=0
#    delta=dot(r1,cosA)
#    
#    return delta
    
if __name__ == '__main__':
    pass
    
