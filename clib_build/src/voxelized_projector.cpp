// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) > (b)) ? (b) : (a))
#define VERY_BIG 1e300
#ifdef WIN32
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <pthread.h>
#include <unistd.h>
#include "DD3_roi_notrans_mm.hpp"
#include <iostream>
#include "getMemorySize.h"

//#define DEBUG
#if defined(DEBUG)

#define DEBUG_00 // Prints "where am I" info
#define DEBUG_10 // Prints info in convert_modular_detector
//#define DEBUG_20 // Prints info in voxelized_projector

#endif

// Global Variable Definitions /////////////

char  TempString[10000];
char  OutputString[10000];
int   PrintReportOutput = 1;

struct module_info
  {
  float *Height;
  float *Width;
  int *Pix;
  float *Coords;
  int *Sub;
  float *Sampling;
  float *Weight;
  float *sourceWeights;
  int nSubSources;
  int maxPixPerModule;
  int maxSubDets;
  int moduleOverlapType;
  int nModuleTypes;
  };

static struct module_info modules = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0,0,0,0,0};

struct material_info
  {
  int materialCount;
  int eBinCount;
  float *muTable;
  };

static struct material_info materials = {0,0,NULL};

static int debug_flag = 0;

struct phantom_info
  {
  int *dims;
  float **volume; // volume is a pointer to pointers because the voxelized phantom may have multiple materials
  float xoff;
  float yoff;
  float zoff;
  float dxy;
  float dz;
  unsigned char *xy_mask;
  };

static struct phantom_info phantom = {NULL,NULL,0,0,0,0,0,NULL};

//////////////////////////////////////////////

void Report(int Unused)
  {
  if (PrintReportOutput)
    {
    sprintf(TempString, "in C> %s", OutputString);
    std::cout << TempString;
    std::cout.flush();
    }
  }

static void dbug(int level,const char *s,...)
  {
  if (debug_flag >= level)
    {
    va_list ap;
    va_start(ap,s);
    vprintf(s,ap);
    va_end(ap);
    fflush(stdout);
    }
}

int* my_memcpyi(int *src, int *dest, int bytes)
  {
  if (dest != NULL) 
    {
    free(dest);
    dest = NULL;
    }
  dest = (int *) malloc(bytes);
  memcpy(dest, src, bytes);
  return dest;
  }

float* my_memcpyf(float *src, float *dest, int bytes)
  {
  if (dest != NULL) 
    {
    free(dest);
    dest = NULL;
    }
  dest = (float *) malloc(bytes);
  memcpy(dest, src, bytes);
  return dest;
  }

extern "C" {
DLLEXPORT void set_src_info_vox(float *sourceWeights, int nSubSources)
  {
  #if defined(DEBUG_00)
  Report(sprintf(OutputString, "In set_src_info_vox\n"));
  #endif

  modules.sourceWeights = my_memcpyf(sourceWeights, modules.sourceWeights, nSubSources*sizeof(float));
  modules.nSubSources = nSubSources;
  }
}

extern "C" {
DLLEXPORT void set_module_info_vox(float *Height, float *Width, int *Pix, float *Coords, int *Sub, float *Sampling, float *Weight, int nModuleTypes, int maxPix, int maxSubDets, int moduleOverlapType)
  {
  int i;
  
  #if defined(DEBUG_00)
  Report(sprintf(OutputString, "In set_module_info_vox\n"));
  #endif

  modules.Height = my_memcpyf(Height, modules.Height, nModuleTypes*sizeof(float));
  modules.Width = my_memcpyf(Width, modules.Width, nModuleTypes*sizeof(float));
  for(i = 0; i < nModuleTypes; i++)
    {
    modules.Height[i] = MAX(1e-7,modules.Height[i]);
    modules.Width[i] = MAX(1e-7,modules.Width[i]);
    }
  modules.Pix = my_memcpyi(Pix, modules.Pix, nModuleTypes*sizeof(int));
  modules.Coords = my_memcpyf(Coords, modules.Coords, maxPix*2*nModuleTypes*sizeof(float));
  modules.Sub = my_memcpyi(Sub, modules.Sub, nModuleTypes*sizeof(int));
  modules.Sampling = my_memcpyf(Sampling, modules.Sampling, 2*maxSubDets*nModuleTypes*sizeof(float));
  modules.Weight = my_memcpyf(Weight, modules.Weight, maxSubDets*nModuleTypes*sizeof(float));
  modules.maxPixPerModule = maxPix;
  modules.maxSubDets = maxSubDets;
  modules.moduleOverlapType = moduleOverlapType;
  modules.nModuleTypes = nModuleTypes;
  }
}

extern "C" {
DLLEXPORT void set_material_info_vox(int materialCount, int eBinCount, float *muTable)
  {
  #if defined(DEBUG_00)
  Report(sprintf(OutputString, "In set_material_info_vox\n"));
  #endif

  materials.materialCount = materialCount;
  materials.eBinCount = eBinCount;
  materials.muTable = my_memcpyf(muTable, materials.muTable, eBinCount*materialCount*sizeof(float));
  }
}

extern "C" {
DLLEXPORT void set_phantom_info_vox(int *Status, float *vol, int *dims, float xoff, float yoff, float zoff, float dxy, float dz, unsigned char *xy_mask, int MaterialIndex, int NumOfMaterials)
  {
  int Material;
  static unsigned long previously_allocated_memory_size;
  unsigned long phantom_num_voxels_one_material  = (unsigned long)dims[0]*(unsigned long)dims[1]*(unsigned long)dims[2];
  unsigned long phantom_num_voxels_all_materials = phantom_num_voxels_one_material*(unsigned long)dims[3];
  unsigned long required_memory_size_one_material  = phantom_num_voxels_one_material *(unsigned long)sizeof(float);
  unsigned long required_memory_size_all_materials = phantom_num_voxels_all_materials*(unsigned long)sizeof(float);
  unsigned long reserved_memory_size = (unsigned long)2*(unsigned long)1024*(unsigned long)1024*(unsigned long)1024; // reserve 2GB memory to store sinogram and other things;
  unsigned long system_memory_size;
  
  *Status = 0;
  
  #if defined(DEBUG_00)
  Report(sprintf(OutputString, "In set_phantom_info_vox\n"));
  #endif

  if (phantom.dims == NULL)
    previously_allocated_memory_size = 0;
  
  if ( (phantom.volume  == NULL) || (previously_allocated_memory_size < required_memory_size_all_materials) )
    {
  	system_memory_size = getMemorySize( );
    Report(sprintf(OutputString, "Preparing to allocate memory for material volume data...\n"));
    Report(sprintf(OutputString, "Phantom is %4i * %4i * %4i voxels * %2i materials,\n", dims[0], dims[1], dims[2], dims[3] ));
    Report(sprintf(OutputString, "         which requires %16lu bytes.\n", required_memory_size_all_materials ));
  	Report(sprintf(OutputString, "  System memory size is %16lu bytes.\n", system_memory_size ));
  	Report(sprintf(OutputString, "              Reserving %16lu bytes for other needs,\n", reserved_memory_size ));
  	Report(sprintf(OutputString, "            this leaves %16lu bytes available for the phantom.\n", (system_memory_size - reserved_memory_size) ));
    if ( (system_memory_size - reserved_memory_size) <  (required_memory_size_all_materials) )
      {
      Report(sprintf(OutputString, "Insuffucient system memory available.\n"));
      *Status = -2;
      return;
      }
    else
      {
      phantom.volume = (float **)malloc(sizeof(float *) * NumOfMaterials);
      if (phantom.volume==NULL)
        {
        Report(sprintf(OutputString, "Memory allocation error - couldn't allocate memory for pointers to materials.\n"));
        *Status = -1;
        return;
        }
      #if defined(DEBUG_00)
      Report(sprintf(OutputString, "Allocated pointers to materials\n"));
      #endif
      for (Material = 1; Material <= NumOfMaterials; Material++ )
        {
        (phantom.volume)[Material-1] = (float *) malloc(required_memory_size_one_material);
        if (phantom.volume[Material-1]==NULL)
          {
          Report(sprintf(OutputString, "Memory allocation error - couldn't allocate memory for material %i.\n", Material));
          *Status = -1;
          return;
          }
        Report(sprintf(OutputString, "Allocated memory for image volume for material %2i\n", Material));
        previously_allocated_memory_size += required_memory_size_one_material;
        }
      Report(sprintf(OutputString, "   Allocated a total of %16lu bytes.\n", previously_allocated_memory_size ));
      }
    }

  Report(sprintf(OutputString, "Copying data for material %2d into C memory...", MaterialIndex));
  memcpy( (phantom.volume)[MaterialIndex-1], vol, dims[0]*dims[1]*dims[2]*sizeof(float));
  Report(sprintf(OutputString,                                                  " done.\n"));
  
  phantom.dims = my_memcpyi(dims, phantom.dims, 4*sizeof(int));
  phantom.xoff = xoff;
  phantom.yoff = yoff;
  phantom.zoff = zoff;
  phantom.dxy = dxy;
  phantom.dz = dz;
  phantom.xy_mask = new unsigned char[dims[0]*dims[1]*2];
  
  memcpy(phantom.xy_mask,xy_mask,sizeof(unsigned char)*dims[0]*dims[1]);
  #if defined(DEBUG_00)
  Report(sprintf(OutputString, "Copied mask into C memory\n"));
  #endif
  unsigned char *transposeMask = phantom.xy_mask + dims[0]*dims[1];
  for(int i = 0;i < dims[0];i++)
    for(int j = 0;j < dims[1];j++)
      transposeMask[j+i*dims[1]] = xy_mask[i+j*dims[0]]; 
  #if defined(DEBUG_00)
  Report(sprintf(OutputString, "Transposed mask\n"));
  #endif
  }
}

int convert_modular_detector(float **xds, float **yds, float **zds, int *nrdetcols, int *nrdetrows, int nModulesIn, int *modTypeInds, float *Up, float *Right, float *Center)
  {
  float tiny = 1e-8, firstU, thisU;
  float tol = 1e-6;
  float *coords;
  int ModTypeIndex, ModIndex, ColIndex, RowIndex;
  int NumPixels, thisNumRows, thisNumCols;
  float *up, *center, *right;
  float inPlane, inZ;

  #if defined(DEBUG_00)
  Report(sprintf(OutputString, "In convert_modular_detector\n"));
  #endif

  //first we verify that each module:
  // 1) has columns oriented in z and rows oriented orth. to z
  // 2) has (u,v) coords that are separable into an xy component and a z component (i.e., a grid)
  // 3) has the same number of rows

  *nrdetrows = 0;
  *nrdetcols = 0;
  // 1) has columns oriented in z and rows oriented orth. to z
  for(int ModIndex = 0; ModIndex < nModulesIn; ModIndex++)
    {
    #if defined(DEBUG_10)
    Report(sprintf(OutputString, "Evaluating module %i of %i for correct orientation\n", ModIndex+1, nModulesIn));
    #endif
    up    = Up    + 3*ModIndex;
    right = Right + 3*ModIndex;
    inPlane = up[0] * up[0] + up[1] * up[1];
    inZ = up[2] * up[2];
    #if defined(DEBUG_10)
    Report(sprintf(OutputString, "up = %f, right = %f, inPlane = %f, inZ = %f\n", *up, *right, inPlane, inZ));
    #endif
    if ((inPlane > inZ) || (inPlane / inZ > tiny) )
      {
      Report(sprintf(OutputString, "ERROR: Modules columns must be parallel to the z direction for the voxelized projector.\n"));
      //dbug(-1,"\r\nERROR: Modules columns must be parallel to the z direction for the voxelized projector.\r\n");
      return(-2);
      }
    inPlane = right[0] * right[0] + right[1] * right[1];
    inZ = right[2] * right[2];
    #if defined(DEBUG_10)
    Report(sprintf(OutputString, "inPlane = %f, inZ = %f\n", inPlane, inZ));
    #endif
    if ( (inZ > inPlane) || (inZ / inPlane > tiny) )
      {
      Report(sprintf(OutputString, "ERROR: Modules rows must be orthogonal to the z direction for the voxelized projector.\n"));
      //dbug(-1,"\r\nERROR: Modules rows must be orthogonal to the z direction for the voxelized projector.\r\n");
      return(-2);
      }
    }
  
  // 2) has (u,v) coords that are separable into an xy component and a z component (i.e., a grid)
  for (ModTypeIndex = 0; ModTypeIndex < modules.nModuleTypes; ModTypeIndex++)
    {
    NumPixels = modules.Pix[ModTypeIndex];
    coords = modules.Coords + ModTypeIndex*modules.maxPixPerModule*2;
    
    #if defined(DEBUG_10)
    Report(sprintf(OutputString, "Evaluating (u,v) coordinates of module type %i of %i\n", ModTypeIndex+1, modules.nModuleTypes));
    Report(sprintf(OutputString, "Module type %i has %i pixels.\n", ModTypeIndex+1, modules.Pix[ModTypeIndex]));
    Report(sprintf(OutputString, "coords[0 1 2 3]: %f %f %f %f \r\n", coords[0], coords[1], coords[2], coords[3]));
    #endif
    //dbug(2,"coords[0 1 2 3]: %f %f %f %f \r\n",coords[0],coords[1],coords[2],coords[3]);
    
    thisNumRows = 1;
    firstU = coords[0];
    thisU  = coords[2*thisNumRows]; // this will be bogus if only one pixel in module but that's ok because
                                    // the loop below will terminate because it's the last pixel.
    // if we're on the same column (same U) and we haven't reached the last pixel, ...
    while ( (fabs(thisU - firstU) < tol) && (thisNumRows != NumPixels) )
      {
      // increment the total number of rows,
      thisNumRows++;
      // and next time through the loop, check the next pixel.
      thisU = coords[2*thisNumRows];
      }

    thisNumCols = modules.Pix[ModTypeIndex] / thisNumRows;
    #if defined(DEBUG_10)
    Report(sprintf(OutputString, "firstU = %f, thisU = %f\n" ,firstU ,thisU));
    #endif
    //dbug(2,"thisU: %f firstU: %f\r\n",thisU,firstU);
    #if defined(DEBUG_10)
    Report(sprintf(OutputString, "thisNumRows = %i, thisNumCols = %i\n", thisNumRows, thisNumCols));
    #endif
    //dbug(2,"pix: %d this_rows: %d\r\n",modules.Pix[ModTypeIndex],thisNumRows);
    
    // 3) has the same number of rows
    if (ModTypeIndex==0)
      *nrdetrows = thisNumRows;
    else if (*nrdetrows != thisNumRows)
      {
      Report(sprintf(OutputString, "ERROR: All module types must have the same number of rows for voxelized projector (%d, %d, %d).\n",
                                                                                                     *nrdetrows, thisNumRows, ModTypeIndex));
      //dbug(-1,"\r\nERROR: All module types must have the same number of rows for voxelized projector (%d, %d, %d).\r\n",*nrdetrows,thisNumRows,ModTypeIndex);
      return(-2);
      }
    //dbug(2,"pix: %d thisNumRows: %d thisNumCols: %d\r\n",modules.Pix[ModTypeIndex],thisNumRows,thisNumCols);
    
    if ((thisNumCols*thisNumRows) != modules.Pix[ModTypeIndex])
      {
      Report(sprintf(OutputString, "ERROR: All module types must be rectilinear grid pixel sampling for voxelized projector.\n"));
      //dbug(-1,"\r\nERROR: All module types must be rectilinear grid pixel sampling for voxelized projector (errorcode 1).\r\n");
      return(-2);
      }

    for (ColIndex = 0; ColIndex < thisNumCols; ColIndex++)
      for (RowIndex = 0; RowIndex < thisNumRows; RowIndex++)
        if ( (fabs(coords[2*(RowIndex+ColIndex*thisNumRows)+1] - coords[2*(RowIndex+0*thisNumRows)+1]) > tol) || 
             (fabs(coords[2*(RowIndex+ColIndex*thisNumRows)]   - coords[2*(0+ColIndex*thisNumRows)]  ) > tol) )
          {
          Report(sprintf(OutputString, "tol: %f \n",tol));
          Report(sprintf(OutputString, "ColIndex: %d RowIndex: %d \r\n",ColIndex,RowIndex));
          Report(sprintf(OutputString, "coords[2*(RowIndex+ColIndex*thisNumRows)+1]: %f coords[2*(RowIndex+0*thisNumRows)+1]: %f \r\n",
                                        coords[2*(RowIndex+ColIndex*thisNumRows)+1],    coords[2*(RowIndex+0*thisNumRows)+1]));
          Report(sprintf(OutputString, "coords[2*(RowIndex+ColIndex*thisNumRows)]:   %f coords[2*(0+ColIndex*thisNumRows)]:   %f \r\n",
                                        coords[2*(RowIndex+ColIndex*thisNumRows)],      coords[2*(0+ColIndex*thisNumRows)]));
          Report(sprintf(OutputString, "ERROR: All module types must be rectilinear grid pixel sampling for voxelized projector.\n"));
          return(-2);
         }
    }

  for (ModIndex = 0; ModIndex < nModulesIn; ModIndex++)
    {
    thisNumCols = modules.Pix[modTypeInds[ModIndex]] / (*nrdetrows);
    *nrdetcols += thisNumCols;    
    }
  #if defined(DEBUG_10)
  Report(sprintf(OutputString, "*nrdetcols = %i, *nrdetrows = %i\n", *nrdetcols, *nrdetrows));
  #endif

  *xds = new float[*nrdetcols];
  *yds = new float[*nrdetcols];
  *zds = new float[*nrdetrows];

  // initial verifications complete... carry on
  // we also verify here that all modules have the same z locations for the rows.
  *nrdetcols = 0;
  for (ModIndex = 0; ModIndex < nModulesIn; ModIndex++)
    {
    #if defined(DEBUG_10)
    Report(sprintf(OutputString, "Assigning system coordinates to pixels of module %i of %i\n", ModIndex+1, nModulesIn));
    #endif
    up     = Up     + 3*ModIndex;
    center = Center + 3*ModIndex;
    right  = Right  + 3*ModIndex;
    #if defined(DEBUG_10)
    Report(sprintf(OutputString, "up = %f, center = %f, right = %f\n", *up, *center, *right));
    #endif
    
    thisNumCols = modules.Pix[modTypeInds[ModIndex]] / (*nrdetrows);
    coords = modules.Coords + modTypeInds[ModIndex]*modules.maxPixPerModule*2;
    if (ModIndex==0)
      {
      for (RowIndex = 0; RowIndex < (*nrdetrows); RowIndex++)
        (*zds)[RowIndex] = center[2] + coords[1+RowIndex*2] * up[2];
      } 
    else
      {
      for (RowIndex = 0; RowIndex < (*nrdetrows); RowIndex++)
        if ( fabs((*zds)[RowIndex] - (center[2] + coords[1+RowIndex*2] * up[2])) > tol )
          {
          Report(sprintf(OutputString, "nzds[RowIndex] = %f, center[2] = %f, coords[1+RowIndex*2] = %f, up[2] = %f, tol = %f\n",
                                      (*zds)[RowIndex],      center[2],      coords[1+RowIndex*2],      up[2],      tol));
          Report(sprintf(OutputString, "ERROR: All modules must line up in z direction (ModIndex = %d, RowIndex = %d).\n", ModIndex, RowIndex));
          return(-2);
          //dbug(0,"\r\nzds[RowIndex]: %f   center[2]: %f   coords[1+RowIndex*2]: %f   up[2]: %f  tol: %f\r\n",(*zds)[RowIndex],center[2],coords[1+RowIndex*2],up[2],tol);
          //dbug(-1,"\r\nERROR: All modules must line up in z direction (ModIndex: %d, RowIndex: %d).\r\n",ModIndex,RowIndex);
          }
      }
    
    for(ColIndex = 0; ColIndex < thisNumCols; ColIndex++)
      {
      (*xds)[ColIndex + *nrdetcols] = center[0] + coords[ColIndex*2*thisNumRows] * right[0];
      (*yds)[ColIndex + *nrdetcols] = center[1] + coords[ColIndex*2*thisNumRows] * right[1];
      }
    
    *nrdetcols += thisNumCols;    
    }
  #if defined(DEBUG_10)
  Report(sprintf(OutputString, "*nrdetcols = %i, *nrdetrows = %i\n", *nrdetcols, *nrdetrows));
  Report(sprintf(OutputString, "Pixel(%3s,%3s) is at (%8s,%8s,%8s):\n", "col", "row", "x", "y", "z"));
  Report(sprintf(OutputString, "Pixel(%3i,%3i) is at (%8.3f,%8.3f,%8.3f)\n",
                                      1,          1,          *xds[1-1],          *yds[1-1],          *zds[1-1]));
  Report(sprintf(OutputString, "Pixel(%3i,%3i) is at (%8.3f,%8.3f,%8.3f)\n",
                                      *nrdetcols, *nrdetrows, *xds[*nrdetcols-1], *yds[*nrdetcols-1], *zds[*nrdetrows-1]));
  #endif

  //dbug(1,"xds[e]: %f\n\r\n",(*xds)[100]);
  
  #if defined(DEBUG_00)
  Report(sprintf(OutputString, "Returning from convert_modular_detector\n"));
  #endif
  return(0);
  }

extern "C" {
DLLEXPORT 
void voxelized_projector(
      int *Status,                    // Output: scalar. 0=normal, 1=Detector definition error
      float unused1,                  // Input: scalar
			float *thisView,                // Output: [nrdetcols*nrdetrows][materials.eBinCount]
			float *sourcePoints,            // Input: [nSubSources][3]
			int nSubSources,                // Input: scalar
			float *unused2,                 // Input: pointer
			int unused3,                    // Input: scalar 
			int *unused4,                   // Input: pointer
			int nModulesIn,                 // Input: scalar
			int *modTypeInds,               // Input: pointer
			float *Up,                      // Input: pointer
			float *Right,                   // Input: pointer
			float *Center,                  // Input: pointer
			int unused5,                    // Input: scalar
			int unused6,                    // Input: scalar
			int MaterialIndex,              // Input: scalar, Material Index to look up mu table 
			int MaterialIndexInMemory,      // Input: scalar, Material Index to materials stored in phantom.volume
			float unused7)                  // Input: scalar
  {
  float *xds = NULL, *yds = NULL, *zds = NULL, x0 = 0, y0 = 0, z0 = 0, totalWt = 0, viewangle;
  int nrdetcols, nrdetrows;

  int EnergyBin;

  *Status = 0; 

  #if defined(DEBUG_00)
  Report(sprintf(OutputString, "In voxelized_projector\n"));
  #endif

  //condense source to one point at center of mass
  for(int i=0;i<nSubSources;i++) {
    x0 += modules.sourceWeights[i] * sourcePoints[0];
    y0 += modules.sourceWeights[i] * sourcePoints[1];
    z0 += modules.sourceWeights[i] * sourcePoints[2];
    totalWt += modules.sourceWeights[i];
    sourcePoints += 3;
  }
  x0 /= totalWt;
  y0 /= totalWt;
  z0 /= totalWt;

  //  viewangle = atan2f(-x0,y0);
  viewangle = 0;

  if (0 != (*Status = convert_modular_detector(&xds, &yds, &zds, &nrdetcols, &nrdetrows, nModulesIn, modTypeInds, Up, Right, Center)))
    {
    Report(sprintf(OutputString, "Error code %i returned by convert_modular_detector\n", *Status));
    return;
    }
      
  float *thisViewProj = new float[nrdetcols*nrdetrows];

  dbug(1,"xds[e]: %f\n\r\n",xds[100]);
  dbug(1,"materials.muTable[0]: %f\r\n",materials.muTable[0]);

  #if defined(DEBUG_00)
  Report(sprintf(OutputString, "Calling DD3Proj_roi_notrans_mm\n"));
  #endif

    float zoffset = 0.0;

  DD3Proj_roi_notrans_mm(x0, y0, z0, nrdetcols, nrdetrows, xds, yds, zds, phantom.xoff, phantom.yoff, phantom.zoff, &viewangle, &zoffset, 1, thisViewProj, phantom.dims[0], phantom.dims[1], phantom.dims[2], (phantom.volume)[MaterialIndexInMemory-1], phantom.dxy, phantom.dz, phantom.xy_mask);

  #if defined(DEBUG_20)
  Report(sprintf(OutputString, "materials.muTable[0] = %f\n", materials.muTable[0]));
  Report(sprintf(OutputString, "thisView[center=%d] = %f\n", nrdetcols*nrdetrows/2,thisView[nrdetcols*nrdetrows/2]));
  #endif

  for(int detIndex=0; detIndex<nrdetcols*nrdetrows; detIndex++)
    for (EnergyBin = 0; EnergyBin < materials.eBinCount; EnergyBin++)
    {
    thisView[detIndex*materials.eBinCount + EnergyBin] = thisViewProj[detIndex]*materials.muTable[EnergyBin * materials.materialCount + (MaterialIndex-1)];
    #if defined(DEBUG_20)
    Report(sprintf(OutputString, "thisView(index1)=%f\n", thisView[detIndex]));
    Report(sprintf(OutputString, "materials.muTable(index2) = %f\n", materials.muTable[EnergyBin * materials.materialCount + (MaterialIndex-1)]));
    #endif
    }
  #if defined(DEBUG_00)
  Report(sprintf(OutputString, "Returning from voxelized_projector\n"));
  #endif

  delete thisViewProj;
  }
}




