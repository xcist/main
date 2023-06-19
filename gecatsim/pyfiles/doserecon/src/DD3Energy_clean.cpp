// GE Proprietary
// last revision: Bruno De Man, Jun 21, 2004

#include <math.h> // fabs(), sqrt(), sin(), cos()
#include <stdlib.h> // malloc(), free()
#include <iostream.h>

extern "C"{
#include "PBDE.hpp"
}

#ifdef WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

#define DETECTOR 0
#define VOLUME   1

//**********************************************************************************************************************
//Selectively comment these out to disable debug output ****************************************************************

// #define DEBUG

#if defined(DEBUG)

char  OutputString[10000];
int   PrintDebug = 1;


#define Report(x) if (PrintDebug) {cout << x; cout.flush();}

#endif

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************

// NOTES:
//  1. All components in the CT system are located by a 3D coordinate system, with the Z-axis at the center of rotation,
//     and the X- and Y-axes in the plane perpendicular to the Z-axis and intersecting the X-ray source
//     when the X-ray source's position is at Z=0.

//  2. All coordinates and distances are in units of mm.

//  3. Parameters related to the X-ray source are named beginning with "Source".

//  4. "Source" parameter names are modified with "X", "Y", and "Z", referring to the above coordinate system.

//  5. Parameters related to the projection domain (projection images are a 2D matrix) are named beginning with
//     "Detector". Some of these parameters relate to system effects which are calculated for each detector pixel,
//     and are named beginning with Detector ase a reminder of the array size and orientation.

//  6. Some "Detector" parameter names are modified with
//        "Col" to refer to columns of detector pixels across the width of the detector (i.e. in the XY plane),
//     or "Row" to refer to rows of detector pixels along the depth of the detector (i.e. in the  Z direction).

//  7. Other "Detector" parameters refer to the absolute position of the detector pixels, within the system.
//     Names of these parameters are modified with "X", "Y", and "Z", referring to the above coordinate system.

//  8. All other parameters relate to the image domain (a 3D volume). These are named beginning with "Volume".

//  9. "Volume" parameter names are modified with "X", "Y", and "Z", referring to the above coordinate system.

// 10. All detector pixels and image volume voxels get projected to the "XZ-plane", i.e. the plane in the image volume
//     where Y = 0.

// 11. Parameters related to the projection (of pixels and voxels) onto the XZ-plane have "Projected" in their names.

// 12. Parameters related to the projection (of pixels and voxels) onto the XZ-plane which pertain to the overlap
//     between detector pixels and image volume voxels are named beginning with "Combined".

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************

#if defined(DEBUG)

void PrintArray(char  *TitleString,
                float *ArrayData_ptr,
                float  DataThreshold,
                int    RowMajor_flag,                     // If true, the X dimension of Array is changing fastest
                int    ArraySize_X,
                int    ArraySize_Y,
                int    IncrementDirection)
  {
  char   PixelCharacter;
  int    Counter_X;
  int    Counter_Y;
  int    Counter_X_Last;
  int    Counter_Y_Last;
  float *Datum_Current_ptr;

  Report(TitleString);
  sprintf(OutputString, "DataThreshold = % 8.3e, RowMajor_flag = % 3i, ArraySize_X = % 3i, ArraySize_Y = % 3i, IncrementDirection = % 3i\n",
                         DataThreshold,          RowMajor_flag,        ArraySize_X,        ArraySize_Y,        IncrementDirection);
  Report(OutputString);
    
  Counter_X_Last = ArraySize_X - 1;
  Counter_Y_Last = ArraySize_Y - 1;
  OutputString[Counter_X_Last + 1] = 0x0A;
  OutputString[Counter_X_Last + 2] = 0x00;

  if (RowMajor_flag)
    {
    Datum_Current_ptr = ArrayData_ptr;
    for (Counter_Y = 0; Counter_Y <= Counter_Y_Last; Counter_Y++)
      {
      for (Counter_X = 0; Counter_X <= Counter_X_Last; Counter_X++)
        {
        if (*Datum_Current_ptr <= DataThreshold)
          PixelCharacter = '.';
        else
          PixelCharacter = '*';
        OutputString[Counter_X] = PixelCharacter;
        Datum_Current_ptr += IncrementDirection;
        }
      Report(OutputString);
      }
    Report("\n");
    }
  else
    {
    for (Counter_Y = 0; Counter_Y <= Counter_Y_Last; Counter_Y++)
      {
      Datum_Current_ptr = ArrayData_ptr + Counter_Y * IncrementDirection;
      for (Counter_X = 0; Counter_X <= Counter_X_Last; Counter_X++)
        {
        if (*Datum_Current_ptr <= DataThreshold)
          PixelCharacter = '.';
        else
          PixelCharacter = '*';
        OutputString[Counter_X] = PixelCharacter;
        Datum_Current_ptr += ArraySize_Y * IncrementDirection;
        }
      Report(OutputString);
      }
    Report("\n");
    }
  }
  
#endif

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
 
// This function populates a vector of pixel boundary coordinates based on a vector of pixel center coordinates.

void DD3Boundaries(int    NumBoundaries,                      // Input:  Scalar
		               float *CenterCoords_ptr,                   // Input:  [NumBoundaries - 1]
		               float *BoundaryCoords_ptr)                 // Output: [NumBoundaries]
  {
  int BoundaryCounter;

  if (NumBoundaries >= 3)
    {
    // The first edge of the first pixel
    *BoundaryCoords_ptr++   = 1.5 * *CenterCoords_ptr - 0.5 * *(CenterCoords_ptr+1);
    for (BoundaryCounter=1 ; BoundaryCounter<=(NumBoundaries-2) ; BoundaryCounter++)
      {
      // Intermediate edges
      *BoundaryCoords_ptr++ = 0.5 * *CenterCoords_ptr + 0.5 * *(CenterCoords_ptr+1);
      CenterCoords_ptr++;
      }
    // The second edge of the last pixel
    *BoundaryCoords_ptr     = 1.5 * *CenterCoords_ptr - 0.5 * *(CenterCoords_ptr-1);
    }
  else
    {
    *(BoundaryCoords_ptr  ) = *CenterCoords_ptr-0.5;
    *(BoundaryCoords_ptr+1) = *CenterCoords_ptr+0.5;
    }
  }
  
//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
 
// This function transposes the first volume, and adds it to the second (untransposed) volume, which is then returned.

void DD3AddTranspose(int    Volume_NumX,                    // Input:  Scalar
		                 int    Volume_NumY,                    // Input:  Scalar
		                 int    Volume_NumZ,                    // Input:  Scalar
		                 float *VolumeToBeTransposed_ptr,       // Input:          [Volume_NumZ * Volume_NumX * Volume_NumY]
		                 float *Volume_NotTransposed_ptr)       // Input & Output: [Volume_NumZ * Volume_NumX * Volume_NumY]
  {
  int Counter_X, Counter_Y, Counter_Z;

  for (Counter_Y=0 ; Counter_Y<=(Volume_NumY-1) ; Counter_Y++)
    {
    for (Counter_X=0 ; Counter_X<=(Volume_NumX-1) ; Counter_X++)
      {
      for (Counter_Z=0 ; Counter_Z<=(Volume_NumZ-1) ; Counter_Z++)
        {
        *Volume_NotTransposed_ptr++ += *VolumeToBeTransposed_ptr++;
        }
      Volume_NotTransposed_ptr += Volume_NumZ * (Volume_NumY-1);
      }
    Volume_NotTransposed_ptr += Volume_NumZ * (1 - Volume_NumX*Volume_NumY);
    }
  }

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************

// Distribute one view of the detector data (as projected to the XZ plane at Y = 0)
// into one row of the image volume (a aingle XZ slab of the volume, projected onto the XZ plane at Y = 0).

void DD3EnergyRow(float  VolumeVoxelBoundaryCoordProjected_X_Initial,  // Input:     Scalar
                  float  VolumeVoxelBoundaryCoordProjected_X_Step,     // Input:     Scalar
                  int    Volume_NumX,                                  // Input:     Scalar
                  float  VolumeVoxelBoundaryCoordProjected_Z_Initial,  // Input:     Scalar
                  float  VolumeVoxelBoundaryCoordProjected_Z_Step,     // Input:     Scalar
                  int    Volume_NumZ,                                  // Input:     Scalar
                  float *Volume_Mu_XZSlab_ptr,                         // Input:     [Volume_NumX][Volume_NumZ]
                  float *Volume_Energy_XZSlab_ptr,                     // Output:    [Volume_NumX][Volume_NumZ]
                  float *DetectorPixelBoundaryCoordsProjected_X_ptr,   // Input:     [Detector_NumCols + 3] 1 sentinal each end
                  int    Detector_IncrementDirection,                  // Input:     Scalar (+1 or -1)
                  float *DetectorPixelBoundaryCoordsShifted_Z_ptr,     // Input:     [Detector_NumRows + 3] 2 sentinals back end
                  float *DetectorMagnificationFactors_ptr,             // Input:     [Detector_NumCols + 2] 1 sentinal each end
                  float  SourceCoord_Z,                                // Input:     Scalar
                  float *Detector_IncidentEnergyViewData_ptr,          // In&Output: [Detector_NumCols    ][Detector_NumRows    ]
                  float *Detector_AbsorbedEnergyViewData_ptr,          // In&Output: [Detector_NumCols    ][Detector_NumRows    ]
                  float *Detector_PathLengths_ptr,                     // Input:     [Detector_NumCols    ][Detector_NumRows    ]
                  int    Detector_NumCols,                             // Input:     Scalar
                  int    Detector_NumRows)                             // Input:     Scalar
  {
  int    FirstActiveBoundary;
  int    JumpCount;
  int    Volume_StayInThisZStack_flag;
  int    DetectorPixelCounter_Col;
  int    DetectorPixelCounter_Col_First;
  int    DetectorPixelCounter_Col_Last;
  int    DetectorPixelCounter_Row;
  int    DetectorPixelCounter_Row_First;
  int    DetectorPixelCounter_Row_Last;
  int    Detector_RowsToMove;
  float  Voxel_XZ_Projected_Area;
  float  Voxel_XZ_Projected_OverlapArea;
  float  Voxel_XZ_Projected_OverlapFraction;
  float  Detector_AbsorbedEnergy_FullVoxel;
  float  Detector_AbsorbedEnergy_OverlapVoxel;
  float  ActiveBoundaryCoordProjected_X;
  float  ActiveBoundaryCoordProjected_Z;
  float  BoundaryCoordProjected_X_Overlap;
  float  BoundaryCoordProjected_Z_Overlap;
  float *DetectorMagnificationFactor_Current_ptr;
  float *DetectorPixelBoundaryCoordProjected_X_Next_ptr;
  float *DetectorPixelBoundaryCoordShifted_Z_Next_ptr;
  float  DetectorPixelBoundaryCoordProjected_Z_Next;
  float *Detector_PathLength_Current_ptr;
  float *Detector_IncidentEnergyViewDatum_Current_ptr;
  float *Detector_AbsorbedEnergyViewDatum_Current_ptr;
  int    VolumeVoxelCounter_X;
  int    VolumeVoxelCounter_X_First;
  int    VolumeVoxelCounter_X_Last;
  int    VolumeVoxelCounter_Z;
  int    VolumeVoxelCounter_Z_First;
  int    VolumeVoxelCounter_Z_Last;
  int    Volume_ZVoxelsToMove;
  float  VolumeVoxelBoundaryCoordProjected_X_Next;
  float  VolumeVoxelBoundaryCoordProjected_Z_Next;
  float *Volume_Mu_CurrentVoxel_ptr;
  float *Volume_Energy_CurrentVoxel_ptr;


  Voxel_XZ_Projected_Area = VolumeVoxelBoundaryCoordProjected_X_Step
                          * VolumeVoxelBoundaryCoordProjected_Z_Step;

  //  Point to the data associated with the first detector pixel and first voxel in the XZ slab.
       DetectorMagnificationFactor_Current_ptr =    DetectorMagnificationFactors_ptr;
  Detector_IncidentEnergyViewDatum_Current_ptr = Detector_IncidentEnergyViewData_ptr;
  Detector_AbsorbedEnergyViewDatum_Current_ptr = Detector_AbsorbedEnergyViewData_ptr;
               Detector_PathLength_Current_ptr =            Detector_PathLengths_ptr;
                    Volume_Mu_CurrentVoxel_ptr =                Volume_Mu_XZSlab_ptr;
                Volume_Energy_CurrentVoxel_ptr =            Volume_Energy_XZSlab_ptr;

  //********************************************************************************************************************
  //  Note that Y is constant; we're dealing with one XZ slab.
  //********************************************************************************************************************
  // Initialize in the X direction:

  if (VolumeVoxelBoundaryCoordProjected_X_Initial > *DetectorPixelBoundaryCoordsProjected_X_ptr)
    {
                               FirstActiveBoundary = VOLUME;
                     VolumeVoxelCounter_X_First    = 0;
           ActiveBoundaryCoordProjected_X          = VolumeVoxelBoundaryCoordProjected_X_Initial;
    DetectorPixelBoundaryCoordProjected_X_Next_ptr = DetectorPixelBoundaryCoordsProjected_X_ptr;
      VolumeVoxelBoundaryCoordProjected_X_Next     = VolumeVoxelBoundaryCoordProjected_X_Initial
                                                   + VolumeVoxelBoundaryCoordProjected_X_Step;
    }
  else
    {
                               FirstActiveBoundary = DETECTOR;
                   DetectorPixelCounter_Col_First  = 0;
           ActiveBoundaryCoordProjected_X          =  *DetectorPixelBoundaryCoordsProjected_X_ptr;
    DetectorPixelBoundaryCoordProjected_X_Next_ptr = ++DetectorPixelBoundaryCoordsProjected_X_ptr;
      VolumeVoxelBoundaryCoordProjected_X_Next     =      VolumeVoxelBoundaryCoordProjected_X_Initial;
      
    }
  
  
  JumpCount = -1; // We are guaranteed at least one jump. If only one, we don't want to advance any of the pointers.
  // Looking along the X direction (parallel to the detector rows),
  if (FirstActiveBoundary == VOLUME)
    {
    // find the first place where a projected detector pixel boundary coordinate
    //            exceeds the first projected volume voxel boundary coordinate.
    // It might be that it takes a number of detector boundaries before that happens
    //   (if the edge of the detector is outside the reconstructed FOV,
    //    or if detector pixels are small and the volume voxels are large).

    while (*DetectorPixelBoundaryCoordProjected_X_Next_ptr <= ActiveBoundaryCoordProjected_X)
      {
      DetectorPixelBoundaryCoordProjected_X_Next_ptr += Detector_IncrementDirection;
      JumpCount++;
      }

    // Synchronize other pointers associated with the detector.
                       DetectorPixelCounter_Col_First = JumpCount;
                                           JumpCount *= Detector_IncrementDirection;
             DetectorMagnificationFactor_Current_ptr += JumpCount;
                                           JumpCount *= Detector_NumRows;
        Detector_IncidentEnergyViewDatum_Current_ptr += JumpCount;
        Detector_AbsorbedEnergyViewDatum_Current_ptr += JumpCount;
                     Detector_PathLength_Current_ptr += JumpCount;
    }
  else
    {
    // find the first place where a projected volume voxel boundary coordinate
    //            exceeds the first projected detector pixel boundary coordinate.
    // It might be that it takes a number of volume boundaries before that happens
    //   (if the edge of the volume is outside the fan beam,
    //    or if voxels are small and the detector pixels are large).

    while (VolumeVoxelBoundaryCoordProjected_X_Next <= ActiveBoundaryCoordProjected_X)
      {
      VolumeVoxelBoundaryCoordProjected_X_Next += VolumeVoxelBoundaryCoordProjected_X_Step;
      JumpCount++;
      }

    // Synchronize other pointers associated with the volume.
        VolumeVoxelCounter_X_First  = JumpCount;
    Volume_Energy_CurrentVoxel_ptr += JumpCount * Volume_NumZ;
        Volume_Mu_CurrentVoxel_ptr += JumpCount * Volume_NumZ;
    }
             
  // X direction initialized.
  //********************************************************************************************************************

  //********************************************************************************************************************
  //  Loop in the X directon over each stack of voxels (stacked in the Z direction).
  
            
          DetectorPixelCounter_Col  = DetectorPixelCounter_Col_First;
              VolumeVoxelCounter_X  =     VolumeVoxelCounter_X_First;
  while ((DetectorPixelCounter_Col <= DetectorPixelCounter_Col_Last)
      && (    VolumeVoxelCounter_X <=     VolumeVoxelCounter_X_Last))
    {
    // The volume voxel counter in the X direction is only conditionally incremented...
    // sometimes this loop is repeated when the detector pixel counter in the X direction is incremented
    // but the volume voxel counter in the X direction is not.
     

    //******************************************************************************************************************
    // Initialize the Z direction:
    
    DetectorPixelBoundaryCoordShifted_Z_Next_ptr      = DetectorPixelBoundaryCoordsShifted_Z_ptr;
    DetectorPixelBoundaryCoordProjected_Z_Next        = SourceCoord_Z 
                                                      + *DetectorPixelBoundaryCoordShifted_Z_Next_ptr++
                                                      * *DetectorMagnificationFactor_Current_ptr;

    if (VolumeVoxelBoundaryCoordProjected_Z_Initial > DetectorPixelBoundaryCoordProjected_Z_Next)
      {
                            FirstActiveBoundary         = VOLUME;
                     VolumeVoxelCounter_Z_First         = 0;
           ActiveBoundaryCoordProjected_Z               = VolumeVoxelBoundaryCoordProjected_Z_Initial;
      VolumeVoxelBoundaryCoordProjected_Z_Next          = VolumeVoxelBoundaryCoordProjected_Z_Initial
                                                        + VolumeVoxelBoundaryCoordProjected_Z_Step;
      }
    else
      {
                                FirstActiveBoundary  = DETECTOR;
                     DetectorPixelCounter_Row_First  = 0;
             ActiveBoundaryCoordProjected_Z          = DetectorPixelBoundaryCoordProjected_Z_Next;
      DetectorPixelBoundaryCoordProjected_Z_Next     = SourceCoord_Z 
                                                     + *DetectorPixelBoundaryCoordShifted_Z_Next_ptr++
                                                     * *DetectorMagnificationFactor_Current_ptr;
        VolumeVoxelBoundaryCoordProjected_Z_Next     = VolumeVoxelBoundaryCoordProjected_Z_Initial;
      
      }

  
    JumpCount = -1; // We are guaranteed at least one jump. If only one, we don't want to advance any of the pointers.
    // Looking along the Z direction (parallel to the detector columns),
    if (FirstActiveBoundary == VOLUME)
      {
      // find the first place where a projected detector pixel boundary coordinate
      //            exceeds the first projected volume voxel boundary coordinate.
      // It might be that it takes a number of detector boundaries before that happens
      //   (if the detector pixels are small and the volume voxels are large).

      while (DetectorPixelBoundaryCoordProjected_Z_Next <= ActiveBoundaryCoordProjected_Z)
        {
        DetectorPixelBoundaryCoordProjected_Z_Next      = SourceCoord_Z 
                                                        + *DetectorPixelBoundaryCoordShifted_Z_Next_ptr++
                                                        * *DetectorMagnificationFactor_Current_ptr;
        JumpCount++;
        }
      
      // Synchronize other variables associated with the detector.

                    DetectorPixelCounter_Row_First  = JumpCount;
                                         JumpCount *= Detector_IncrementDirection;
      Detector_IncidentEnergyViewDatum_Current_ptr += JumpCount;
      Detector_AbsorbedEnergyViewDatum_Current_ptr += JumpCount;
                   Detector_PathLength_Current_ptr += JumpCount;
      }
    else
      {
      // find the first place where a projected volume voxel boundary coordinate
      //            exceeds the first projected detector pixel boundary coordinate.
      // It might be that it takes a number of volume boundaries before that happens
      //   (if the edge of the volume is outside the cone beam,
      //    or if voxels are small and the detector pixels are large).

      while (VolumeVoxelBoundaryCoordProjected_Z_Next <= ActiveBoundaryCoordProjected_Z)
        {
        VolumeVoxelBoundaryCoordProjected_Z_Next += VolumeVoxelBoundaryCoordProjected_Z_Step;
        JumpCount++;
        }

      // Synchronize other variables associated with the volume.
          VolumeVoxelCounter_Z_First  = JumpCount;
      Volume_Energy_CurrentVoxel_ptr += JumpCount;
          Volume_Mu_CurrentVoxel_ptr += JumpCount;
      }

    // Z direction initialized.
    //******************************************************************************************************************
    
    //******************************************************************************************************************
    // Looking along the X direction (parallel to the detector rows),
    // find the amount of X-overlap between the boundary currently under consideration, and the next boundary.

    if (VolumeVoxelBoundaryCoordProjected_X_Next <= *DetectorPixelBoundaryCoordProjected_X_Next_ptr)
      {
      // Next X boundary is a volume boundary. Calculate the X overlap.
        
                    BoundaryCoordProjected_X_Overlap   = VolumeVoxelBoundaryCoordProjected_X_Next
                                                            - ActiveBoundaryCoordProjected_X;
                         
      // Next time through, we'll move to the next Z-stack in the XZ slab.
         
      Volume_StayInThisZStack_flag = 0;
      
      // Increment to the next volume voxel (in the X direction).
                         
                        VolumeVoxelCounter_X++;
              ActiveBoundaryCoordProjected_X           = VolumeVoxelBoundaryCoordProjected_X_Next;
         VolumeVoxelBoundaryCoordProjected_X_Next     += VolumeVoxelBoundaryCoordProjected_X_Step;
         
      }
    else
      {
      // Next X boundary is a detector boundary. Calculate the X overlap.
        
                    BoundaryCoordProjected_X_Overlap   = *DetectorPixelBoundaryCoordProjected_X_Next_ptr
                                                               - ActiveBoundaryCoordProjected_X;
                    
      // Next time through, we'll stay in this Z-stack in the XZ slab, to go through it again with the
      // next column of detector pixels, at least part of which overlaps with with this Z-stack in the XZ slab.

      Volume_StayInThisZStack_flag = 1;
      
      // Increment to the next detector pixel (in the X direction).
                         
                      DetectorPixelCounter_Col++;
              ActiveBoundaryCoordProjected_X           = *DetectorPixelBoundaryCoordProjected_X_Next_ptr;
       DetectorPixelBoundaryCoordProjected_X_Next_ptr += Detector_IncrementDirection;
              DetectorMagnificationFactor_Current_ptr += Detector_IncrementDirection;
              
    }
    
    // X-overlap found.
    //******************************************************************************************************************
    
    //******************************************************************************************************************
    // Note that both Y and X are now constant; we're dealing with one Z-stack (voxels stacked in the Z direction).
    // Looking along the Z direction (parallel to the detector columns),
    // find the amount of Z-overlap between the boundary currently under consideration, and the next boundary.

    // Then, based on the Z-overlap and the X-overlap (previously found), project a fraction of the current detector
    // pixel data into the current image volume voxel.
        
              VolumeVoxelCounter_Z    =     VolumeVoxelCounter_Z_First;
            DetectorPixelCounter_Row  = DetectorPixelCounter_Row_First;
            

    while ((  VolumeVoxelCounter_Z   <=     VolumeVoxelCounter_Z_Last)
        && (DetectorPixelCounter_Row <= DetectorPixelCounter_Row_Last))
      {
      // The volume voxel counter in the Z direction is only conditionally incremented...
      // sometimes this loop is repeated when the detector pixel counter in the Z direction is incremented
      // but the volume voxel counter in the Z direction is not.
    
      if (VolumeVoxelBoundaryCoordProjected_Z_Next <= DetectorPixelBoundaryCoordProjected_Z_Next) 
        {
        // Next Z boundary is a volume boundary. Calculate the Z overlap.
          
        BoundaryCoordProjected_Z_Overlap  = VolumeVoxelBoundaryCoordProjected_Z_Next
                                          - ActiveBoundaryCoordProjected_Z;
        }
      else 
        {
        // Next Z boundary is a detector boundary. Calculate the Z overlap,
          
        BoundaryCoordProjected_Z_Overlap = DetectorPixelBoundaryCoordProjected_Z_Next
                                         - ActiveBoundaryCoordProjected_Z;
        }

      // Calculate the overlap area,
      Voxel_XZ_Projected_OverlapArea = BoundaryCoordProjected_X_Overlap
                                     * BoundaryCoordProjected_Z_Overlap;

      Voxel_XZ_Projected_OverlapFraction = Voxel_XZ_Projected_OverlapArea
                                         / Voxel_XZ_Projected_Area;

      // the energy associated with the full voxel area,
      Detector_AbsorbedEnergy_FullVoxel = *Detector_IncidentEnergyViewDatum_Current_ptr * (1 - exp(-*Volume_Mu_CurrentVoxel_ptr * *Detector_PathLength_Current_ptr));
      
      // and the energy associated with the overlap area.
      Detector_AbsorbedEnergy_OverlapVoxel = Detector_AbsorbedEnergy_FullVoxel * Voxel_XZ_Projected_OverlapFraction;
      
      // Add that energy to the absorbed energy volume, as well as the detector array
      //                                               (to be subtracted from the incident energy detector array later).
                    *Volume_Energy_CurrentVoxel_ptr += Detector_AbsorbedEnergy_OverlapVoxel;
      *Detector_AbsorbedEnergyViewDatum_Current_ptr += Detector_AbsorbedEnergy_OverlapVoxel;
      

      if (VolumeVoxelBoundaryCoordProjected_Z_Next <= DetectorPixelBoundaryCoordProjected_Z_Next) 
        {
        // Next Z boundary is a volume boundary. Increment to the next volume voxel (in the Z direction).
                                          
                       VolumeVoxelCounter_Z++;
             ActiveBoundaryCoordProjected_Z       = VolumeVoxelBoundaryCoordProjected_Z_Next;
        VolumeVoxelBoundaryCoordProjected_Z_Next += VolumeVoxelBoundaryCoordProjected_Z_Step;
            Volume_Mu_CurrentVoxel_ptr++;
        Volume_Energy_CurrentVoxel_ptr++;
        
#if defined(DEBUG_29) || defined(DEBUG_39) || defined(DEBUG_40)
        if (PrintZInfo)
          {
          sprintf(OutputString, "                                              Moving to the next Z voxel in the    stack: DetectorPixelCounter_Row = %3i, VolumeVoxelCounter_Z = %3i\n", DetectorPixelCounter_Row, VolumeVoxelCounter_Z);
          Report(OutputString);
          }
#endif
        }
      else 
        {
        // Next Z boundary is a detector boundary. Increment to the next detector pixel (in the Z direction).
                                          
                          DetectorPixelCounter_Row++;
                  ActiveBoundaryCoordProjected_Z      = DetectorPixelBoundaryCoordProjected_Z_Next;
           DetectorPixelBoundaryCoordProjected_Z_Next = SourceCoord_Z
                                                      + *(++DetectorPixelBoundaryCoordShifted_Z_Next_ptr)
                                                      * *DetectorMagnificationFactor_Current_ptr;
        Detector_IncidentEnergyViewDatum_Current_ptr += Detector_IncrementDirection;
        Detector_AbsorbedEnergyViewDatum_Current_ptr += Detector_IncrementDirection;
                     Detector_PathLength_Current_ptr += Detector_IncrementDirection;
                     
#if defined(DEBUG_29) || defined(DEBUG_39) || defined(DEBUG_40)
        if (PrintZInfo)
          {
          sprintf(OutputString, "                                              Moving to the next Z pixel in the detector: DetectorPixelCounter_Row = %3i, VolumeVoxelCounter_Z = %3i\n", DetectorPixelCounter_Row, VolumeVoxelCounter_Z);
          Report(OutputString);
          }
#endif
        }
      } // End Z-loop.
      //****************************************************************************************************************
    
    if (Volume_StayInThisZStack_flag)
      {
      // Continue on to the next detector column
                               Detector_RowsToMove  = (Detector_NumRows - DetectorPixelCounter_Row) * Detector_IncrementDirection;
      Detector_IncidentEnergyViewDatum_Current_ptr += Detector_RowsToMove;
      Detector_AbsorbedEnergyViewDatum_Current_ptr += Detector_RowsToMove;
                   Detector_PathLength_Current_ptr += Detector_RowsToMove;
                   
      // but go back to the beginning of this Z-stack in the XZ slab.
                              Volume_ZVoxelsToMove  = VolumeVoxelCounter_Z;
                    Volume_Energy_CurrentVoxel_ptr -= Volume_ZVoxelsToMove;
                        Volume_Mu_CurrentVoxel_ptr -= Volume_ZVoxelsToMove;
                        
#if defined(DEBUG_29) || defined(DEBUG_35) || defined(DEBUG_39) || defined(DEBUG_40)
      sprintf(OutputString, "Moving to the next            detector column. Moving forward  %3i detector rows.\n", Detector_RowsToMove * Detector_IncrementDirection);
      Report(OutputString);
      sprintf(OutputString, "Staying in the current Z-stack in the XZ slab. Moving backward %3i voxels.\n",        Volume_ZVoxelsToMove);
      Report(OutputString);
#endif
      }
    else
      {
      // Continue on to the next Z-stack in the XZ slab,
                              Volume_ZVoxelsToMove  = Volume_NumZ - VolumeVoxelCounter_Z;
                    Volume_Energy_CurrentVoxel_ptr += Volume_ZVoxelsToMove;
                        Volume_Mu_CurrentVoxel_ptr += Volume_ZVoxelsToMove;
          
      // but go back to the beginning of this detector column.
                               Detector_RowsToMove  = DetectorPixelCounter_Row*Detector_IncrementDirection;
      Detector_IncidentEnergyViewDatum_Current_ptr -= Detector_RowsToMove;
      Detector_AbsorbedEnergyViewDatum_Current_ptr -= Detector_RowsToMove;
                   Detector_PathLength_Current_ptr -= Detector_RowsToMove;
                   
#if defined(DEBUG_29) || defined(DEBUG_35) || defined(DEBUG_39) || defined(DEBUG_40)
      sprintf(OutputString, "Moving to the next     Z-stack in the XZ slab. Moving forward  %3i voxels.\n",        Volume_ZVoxelsToMove);
      Report(OutputString);
      sprintf(OutputString, "Staying in the current        detector column. Moving backward %3i detector rows.\n", Detector_RowsToMove * Detector_IncrementDirection);
      Report(OutputString);
#endif
      }
    
    } // End X-loop.
    //******************************************************************************************************************
 
  // Subtract the energy absorbed by this XZ slab from the what remains of the incident energy.
  
#if defined(DEBUG_29) || defined(DEBUG_40)
  Report("Subtract the energy absorbed by this XZ slab from the what remains of the incident energy.\n");
#endif

  int    DetectorPixelCounter;
  int    DetectorPixelCounter_Last;
  
  DetectorPixelCounter_Last = Detector_NumCols * Detector_NumRows - 1;
  
  Detector_IncidentEnergyViewDatum_Current_ptr = Detector_IncidentEnergyViewData_ptr;
  Detector_AbsorbedEnergyViewDatum_Current_ptr = Detector_AbsorbedEnergyViewData_ptr;
  
  for (DetectorPixelCounter = 0; DetectorPixelCounter <= DetectorPixelCounter_Last; DetectorPixelCounter++)
    {
    *Detector_IncidentEnergyViewDatum_Current_ptr -= *Detector_AbsorbedEnergyViewDatum_Current_ptr;
     Detector_IncidentEnergyViewDatum_Current_ptr += Detector_IncrementDirection;
     Detector_AbsorbedEnergyViewDatum_Current_ptr += Detector_IncrementDirection;
    }
  }

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
 
void DD3EnergyView(float  SourceCoord_X,                              // Input:  Scalar
                   float  SourceCoord_Y,                              // Input:  Scalar
                   float  SourceCoord_Z,                              // Input:  Scalar
                   int    Detector_NumCols,                           // Input:  Scalar
                   int    Detector_NumRows,                           // Input:  Scalar
                   int    ProjectionVertical_flag,                    // Input:  Scalar
                   float *DetectorPixelBoundaryCoordsRotated_X_ptr,   // Input:  [Detector_NumCols + 1]
                   float *DetectorPixelBoundaryCoordsRotated_Y_ptr,   // Input:  [Detector_NumRows + 1]
                   float *DetectorPixelBoundaryCoordsProjected_X_ptr, // Output: [Detector_NumCols + 3] 1 sentinal each end
                   float *DetectorPixelBoundaryCoordsShifted_Z_ptr,   // Input:  [Detector_NumRows + 3] 2 sentinals back end
                   float *DetectorMagnificationFactors_ptr,           // Output: [Detector_NumCols + 2] 1 sentinal each end
                   float  DetectorPixelAspectRatio_Z_X,               // Input:  Scalar
                   float *Detector_IncidentEnergyViewData_ptr,        // Input:  [Detector_NumCols  ][Detector_NumRows  ]
                   float *Detector_AbsorbedEnergyViewData_ptr,        // Output: [Detector_NumCols  ][Detector_NumRows  ]
                   float *Detector_PathLengths_ptr,                   // Output: [Detector_NumCols  ][Detector_NumRows  ]
                   int    Volume_NumX,                                // Input:  Scalar
                   int    Volume_NumY,                                // Input:  Scalar
                   int    Volume_NumZ,                                // Input:  Scalar
                   float  Volume_VoxelSize_X,                         // Input:  Scalar
                   float  Volume_VoxelSize_Y,                         // Input:  Scalar
                   float  Volume_VoxelSize_Z,                         // Input:  Scalar
                   float *Volume_Mu_Untransposed_ptr,                 // Input:  [Volume_NumY][Volume_NumX][Volume_NumZ]
                   float *Volume_Mu_Transposed_ptr,                   // Input:  [Volume_NumY][Volume_NumX][Volume_NumZ]
                   float *Volume_Energy_Untransposed_ptr,             // Output: [Volume_NumY][Volume_NumX][Volume_NumZ]
                   float *Volume_Energy_Transposed_ptr)               // Output: [Volume_NumY][Volume_NumX][Volume_NumZ]
  {
  int    Temp_int;
  int    Counter;
  float  Temp_float;
  float  CosineFactor;
  int    Detector_NumPixels;
  int    Detector_NumPixelBoundaries_XY;
  int    DetectorCounter_Row;
  int    DetectorCounter_Row_Last;
  int    DetectorCounter_Col;
  int    DetectorCounter_Col_Last;
  int    DetectorPixelBoundaryCounter_XY_Last;
  int    Detector_IncrementDirection;
  float *DetectorMagnificationFactor_Current_ptr;
  float *DetectorPixelBoundaryCoordRotated_X_Current_ptr;
  float *DetectorPixelBoundaryCoordRotated_Y_Current_ptr;
  float *DetectorPixelBoundaryCoordsProjected_X_Initial_ptr;
  float *DetectorPixelBoundaryCoordProjected_X_Current_ptr;
  float *DetectorPixelBoundaryCoordShifted_Z_Current_ptr;
  float  DetectorPixelBoundaryCoordsProjected_X_Midpoint;
  float  DetectorPixelBoundaryCoordsProjected_X_Step;
  float  DetectorPixelBoundaryCoordsShifted_Z_Midpoint;
  float  DetectorPixelBoundaryCoordsShifted_Z_Step;
  float *Detector_PathLength_Current_ptr;
  float *Detector_IncidentEnergyViewData_Current_ptr;
  float *Detector_AbsorbedEnergyViewData_Current_ptr;
  int    Volume_NumVoxels_XZSlab;
  int    VolumeCounter_Y;
  int    VolumeCounter_Y_First;
  int    VolumeCounter_Y_Last;
  int    VolumeCounter_Y_BeyondLast;
  int    VolumeCounter_Y_IncrementDirection;
  float  VolumeVoxelBoundaryCoordProjected_X_Initial;
  float  VolumeVoxelBoundaryCoordProjected_X_Step;
  float  VolumeVoxelBoundaryCoordProjected_Z_Initial;
  float  VolumeVoxelBoundaryCoordProjected_Z_Step;
  float  VolumeXZSlab_Y_Current;
  float  VolumeXZSlabMagnificationFactor_Current;
  float *Volume_Mu_Active_ptr;
  float *Volume_Mu_XZSlab_Current_ptr;
  float *Volume_Energy_Active_ptr;
  float *Volume_Energy_XZSlab_Current_ptr;

  Detector_NumPixels                   = Detector_NumCols * Detector_NumRows;
    Volume_NumVoxels_XZSlab            =   Volume_NumX    *   Volume_NumZ;

  // C array referencing 0 -> n-1
  DetectorCounter_Col_Last             = Detector_NumCols - 1;
  DetectorCounter_Row_Last             = Detector_NumRows - 1;
  
  // 1 more boundary than pixels
  Detector_NumPixelBoundaries_XY       = Detector_NumCols + 1;
  // C array referencing 0 -> n-1
  DetectorPixelBoundaryCounter_XY_Last = Detector_NumPixelBoundaries_XY - 1;
  
#if defined(DEBUG_50)
  PrintArray("\nIncident energy projection:\n",   Detector_IncidentEnergyViewData_ptr, 1e-4, 0, Detector_NumCols, Detector_NumRows, 1);

  PrintArray("\nAbsorbed energy projection:\n",   Detector_AbsorbedEnergyViewData_ptr, 0.0,  0, Detector_NumCols, Detector_NumRows, 1);

  for (VolumeCounter_Y = 0; VolumeCounter_Y < Volume_NumY/2; VolumeCounter_Y++)
    {
    sprintf(OutputString, "\nMu Volume XZ slab %3i:\n",      VolumeCounter_Y);
    PrintArray(OutputString, Volume_Mu_Untransposed_ptr    + VolumeCounter_Y * Volume_NumVoxels_XZSlab,           0.0,  0, Volume_NumX,      Volume_NumZ,      1);
    }

  for (VolumeCounter_Y = 0; VolumeCounter_Y < Volume_NumY/2; VolumeCounter_Y++)
    {
    sprintf(OutputString, "\nEnergy Volume XZ slab %3i:\n",  VolumeCounter_Y);
    PrintArray(OutputString, Volume_Energy_Untransposed_ptr + VolumeCounter_Y * Volume_NumVoxels_XZSlab,           0.0,  0, Volume_NumX,      Volume_NumZ,      1);
    }
#endif

  // *******************************************************************************************************************
  // Project the X coordinates of the detector pixel boundaries onto the X axis.
  
#if defined(DEBUG_00)
  Report("Project the X coordinates of the detector pixel boundaries onto the X axis.\n");
#endif
  // Initialize:
  // Pointers to the first X-boundary coordinate, 
  DetectorPixelBoundaryCoordRotated_X_Current_ptr   = DetectorPixelBoundaryCoordsRotated_X_ptr;
  DetectorPixelBoundaryCoordRotated_Y_Current_ptr   = DetectorPixelBoundaryCoordsRotated_Y_ptr;
  // and pointer to the associated magnification factor.
  DetectorMagnificationFactor_Current_ptr = DetectorMagnificationFactors_ptr + 1; // +1 to skip the sentinal.
  
  DetectorPixelBoundaryCoordProjected_X_Current_ptr = DetectorPixelBoundaryCoordsProjected_X_ptr + 1;  // +1 to skip the sentinal.  
  if (ProjectionVertical_flag)
    {
    // We'll be doing a vertical backprojection.
      
    for (Counter=0 ; Counter <= DetectorPixelBoundaryCounter_XY_Last; Counter++)
      {
      // Calculate the X coordinate.
        
      *DetectorPixelBoundaryCoordProjected_X_Current_ptr++ =  (SourceCoord_X * *DetectorPixelBoundaryCoordRotated_Y_Current_ptr
                                                             - SourceCoord_Y * *DetectorPixelBoundaryCoordRotated_X_Current_ptr)
                                                            / (*DetectorPixelBoundaryCoordRotated_Y_Current_ptr - SourceCoord_Y);
      *DetectorMagnificationFactor_Current_ptr++ = SourceCoord_Y
                                                 /(SourceCoord_Y - *DetectorPixelBoundaryCoordRotated_Y_Current_ptr);
      DetectorPixelBoundaryCoordRotated_Y_Current_ptr++;
      DetectorPixelBoundaryCoordRotated_X_Current_ptr++;
      }
    
    // We will work with the "untransposed" volume arrays.
      
    Volume_Mu_Active_ptr     = Volume_Mu_Untransposed_ptr;
    Volume_Energy_Active_ptr = Volume_Energy_Untransposed_ptr;
    }
  else
    {
    // We'll be doing a horizontal backprojection.
      
    for (Counter=0 ; Counter <= DetectorPixelBoundaryCounter_XY_Last ; Counter++)
      {
      // Calculate the y-intercept in the row direction: i.e. negative y-axis, so after transpose = positive x-axis.
      
      *DetectorPixelBoundaryCoordProjected_X_Current_ptr++ = -(SourceCoord_Y * *DetectorPixelBoundaryCoordRotated_X_Current_ptr
                                                             - SourceCoord_X * *DetectorPixelBoundaryCoordRotated_Y_Current_ptr)
                                                            / (*DetectorPixelBoundaryCoordRotated_X_Current_ptr - SourceCoord_X);
      // Calculate the magnification factor.
      
      *DetectorMagnificationFactor_Current_ptr++ = SourceCoord_X 
                                                 /(SourceCoord_X - *DetectorPixelBoundaryCoordRotated_X_Current_ptr);
      DetectorPixelBoundaryCoordRotated_X_Current_ptr++;
      DetectorPixelBoundaryCoordRotated_Y_Current_ptr++;
      }
    
      // We will work with the "transposed" volume arrays.
    
    Volume_Mu_Active_ptr     = Volume_Mu_Transposed_ptr;
    Volume_Energy_Active_ptr = Volume_Energy_Transposed_ptr;
    
    // We also need to transpose a few related variables.
    
    Temp_int      = Volume_NumY;
    Volume_NumY   = Volume_NumX;
    Volume_NumX   = Temp_int;
    Temp_float    =  SourceCoord_Y;
    SourceCoord_Y = -SourceCoord_X;
    SourceCoord_X = -Temp_float;
    }
 
  // Done projecting the X coordinates of the detector pixel boundaries onto the X axis.
  // *******************************************************************************************************************
  // *******************************************************************************************************************
  // Seems like this is adjusting the magnification factors values by averaging each with the next.
  // Not sure why this needs to happen - because the calculation of these factors was based on edges, not centers?

  // Copy the 2nd element (the 1st "real" value) into the 1st element (the 1st sentinel).
  *DetectorMagnificationFactors_ptr = *(DetectorMagnificationFactors_ptr + 1);
  // Starting with the 2nd element (the 1st "real" value),
  //   ending with the 2nd-to-last element (the last "real" value),
  // average each value with the next.
  for (DetectorMagnificationFactor_Current_ptr =  DetectorMagnificationFactors_ptr + 1;
       DetectorMagnificationFactor_Current_ptr <= DetectorMagnificationFactors_ptr + Detector_NumCols;
       DetectorMagnificationFactor_Current_ptr++)
    {
    *DetectorMagnificationFactor_Current_ptr = (*(DetectorMagnificationFactor_Current_ptr    )
                                              + *(DetectorMagnificationFactor_Current_ptr + 1))
                                              /2.;
    }

  // Done averaging the magnification factors.
  // *******************************************************************************************************************
  // *******************************************************************************************************************
  // Calculate the path length through a voxel slab, on a line from from the source to each detector pixel.
  // These are independent of volume XZ slab.
  
#if defined(DEBUG_00)
  Report("Calculate the path length through a voxel slab, on a line from from the source to each detector pixel.\n");
  Report("These are independent of volume XZ slab.\n");
#endif

  // Initialize:
  // Pointer to the first path length,

  Detector_PathLength_Current_ptr = Detector_PathLengths_ptr;

  // pointer to the X coordinate of the first boundary of the first "real" pixel in the detector,

  DetectorPixelBoundaryCoordProjected_X_Current_ptr = DetectorPixelBoundaryCoordsProjected_X_ptr + 1;  // +1 for the sentinal.
  
  // and pointer to the associated magnification factor.
  
  DetectorMagnificationFactor_Current_ptr = DetectorMagnificationFactors_ptr + 1; // +1 for the sentinal.
  
  for (DetectorCounter_Col=0 ; DetectorCounter_Col<=DetectorCounter_Col_Last ; DetectorCounter_Col++)
    {
    // Some things are the same for all rows in the column.
    
    // Find the X coordinate of the midpoint between the current and the next column.
      
    DetectorPixelBoundaryCoordsProjected_X_Midpoint = (*(DetectorPixelBoundaryCoordProjected_X_Current_ptr  )
                                                     + *(DetectorPixelBoundaryCoordProjected_X_Current_ptr+1))
                                                        /2.
                                                    - SourceCoord_X;
    
    // Find the X distance between the current and the next column.
    
    DetectorPixelBoundaryCoordsProjected_X_Step = fabs(*(DetectorPixelBoundaryCoordProjected_X_Current_ptr+1)
                                                     - *(DetectorPixelBoundaryCoordProjected_X_Current_ptr  ));
    
    // The Z coordinates are the same for all columns but vary by row. Initialize to the first row.
    
    DetectorPixelBoundaryCoordShifted_Z_Current_ptr = DetectorPixelBoundaryCoordsShifted_Z_ptr;
    
    for (DetectorCounter_Row=0 ; DetectorCounter_Row<=DetectorCounter_Row_Last ; DetectorCounter_Row++)
      {
      // Find the Z coordinate of the midpoint between the current row boundary and the next row boundary.

      DetectorPixelBoundaryCoordsShifted_Z_Midpoint = *DetectorMagnificationFactor_Current_ptr
                                                      * (*(DetectorPixelBoundaryCoordShifted_Z_Current_ptr  )
                                                       + *(DetectorPixelBoundaryCoordShifted_Z_Current_ptr+1))
                                                          /2.;

      // Find the Z distance between the current row boundary and the next row boundary.

      DetectorPixelBoundaryCoordsShifted_Z_Step     = *DetectorMagnificationFactor_Current_ptr
                                                  * fabs(*(DetectorPixelBoundaryCoordShifted_Z_Current_ptr+1)
                                                       - *(DetectorPixelBoundaryCoordShifted_Z_Current_ptr  ));
      // Calculate the path length.

      CosineFactor = fabs(SourceCoord_Y)                                   // adjacent
                                                                           // hypotenuse by 3D pythagoreum:
                   / sqrt((SourceCoord_Y                                   // adjacent^2
                         * SourceCoord_Y)
                        + (DetectorPixelBoundaryCoordsProjected_X_Midpoint // opposite^2
                         * DetectorPixelBoundaryCoordsProjected_X_Midpoint)
                                                                           // the other opposite^2
                        + (DetectorPixelBoundaryCoordsShifted_Z_Midpoint * DetectorPixelAspectRatio_Z_X
                         * DetectorPixelBoundaryCoordsShifted_Z_Midpoint * DetectorPixelAspectRatio_Z_X));

      *Detector_PathLength_Current_ptr++ = Volume_VoxelSize_Y / CosineFactor;

      DetectorPixelBoundaryCoordShifted_Z_Current_ptr++;
        
#if defined(DEBUG_05)
      if ((DetectorCounter_Col == DetectorCounter_Col_Last/2) && (DetectorCounter_Row == DetectorCounter_Row_Last/2))
        {
        sprintf(OutputString, "DetectorCounter_Col = %3i, DetectorCounter_Row = %3i, Detector_PathLength_Current = % 8.2f\n",
                               DetectorCounter_Col,       DetectorCounter_Row,     *(Detector_PathLength_Current_ptr-1));
        Report(OutputString);
        }
#endif
      } // end detector row loop
    
    DetectorPixelBoundaryCoordProjected_X_Current_ptr++;
    DetectorMagnificationFactor_Current_ptr++;
    
    } // end detector column loop

  // Set the sentinels.

    *(DetectorPixelBoundaryCoordsProjected_X_ptr                         )  // First sentinel
  = *(DetectorPixelBoundaryCoordsProjected_X_ptr                      + 1); // = first "real" value.

    *(DetectorPixelBoundaryCoordsProjected_X_ptr + Detector_NumCols+1 + 1)  // Last sentinel
  = *(DetectorPixelBoundaryCoordsProjected_X_ptr + Detector_NumCols+1    ); // = last "real" value.

  // Done projecting view onto the XZ plane.
  // *******************************************************************************************************************
  // *******************************************************************************************************************
  // Initialize direction to progress in the X direction.

#if defined(DEBUG_00)
  Report("Initialize direction to progress in the X direction.\n");
#endif

  if (*(DetectorPixelBoundaryCoordsProjected_X_ptr+2)
    > *(DetectorPixelBoundaryCoordsProjected_X_ptr+1))
    {
    // Go forward.
    
    Detector_IncrementDirection = 1;
    
    // Point to the first element of the arrays.
    
#if defined(DEBUG_19)
       DetectorPixelBoundaryCoordRotated_X_Current_ptr = DetectorPixelBoundaryCoordsRotated_X_ptr;
       DetectorPixelBoundaryCoordRotated_Y_Current_ptr = DetectorPixelBoundaryCoordsRotated_Y_ptr;
#endif
    DetectorPixelBoundaryCoordsProjected_X_Initial_ptr = DetectorPixelBoundaryCoordsProjected_X_ptr + 1; // skip the sentinel
               DetectorMagnificationFactor_Current_ptr =           DetectorMagnificationFactors_ptr;
                       Detector_PathLength_Current_ptr =                   Detector_PathLengths_ptr;
           Detector_IncidentEnergyViewData_Current_ptr =        Detector_IncidentEnergyViewData_ptr;
           Detector_AbsorbedEnergyViewData_Current_ptr =        Detector_AbsorbedEnergyViewData_ptr;
    }  
  else
    {
    // Go backward.
    
    Detector_IncrementDirection = -1;
    
    // Point to the last element of the arrays.
    //                                                                   POINTER ARITHMETIC NOTES: +1 for one more boundary than pixels,
    //                                                                                             +1 or +2 for sentinals,
    //                                                                                             -1 to get back to last element,
    //                                                                                             -2 to skip the last sentinal (get back to last "real" value)
#if defined(DEBUG_19)
    DetectorPixelBoundaryCoordRotated_X_Current_ptr    = DetectorPixelBoundaryCoordsRotated_X_ptr   + (Detector_NumCols+1    )                        - 1;
    DetectorPixelBoundaryCoordRotated_Y_Current_ptr    = DetectorPixelBoundaryCoordsRotated_Y_ptr   + (Detector_NumCols+1    )                        - 1;
#endif
    DetectorPixelBoundaryCoordsProjected_X_Initial_ptr = DetectorPixelBoundaryCoordsProjected_X_ptr + (Detector_NumCols+1 + 2)                        - 2;
               DetectorMagnificationFactor_Current_ptr =           DetectorMagnificationFactors_ptr + (Detector_NumCols   + 2)                        - 1;
                       Detector_PathLength_Current_ptr =                   Detector_PathLengths_ptr + (Detector_NumCols      )*(Detector_NumRows    ) - 1;
           Detector_IncidentEnergyViewData_Current_ptr =        Detector_IncidentEnergyViewData_ptr + (Detector_NumCols      )*(Detector_NumRows    ) - 1;
           Detector_AbsorbedEnergyViewData_Current_ptr =        Detector_AbsorbedEnergyViewData_ptr + (Detector_NumCols      )*(Detector_NumRows    ) - 1;
    }
  // Done initializing direction to progress in the X direction.
  // *******************************************************************************************************************
  // *******************************************************************************************************************
  // Initialize direction to progress in the Y direction.
  
  if (SourceCoord_Y > 0)
    {
    VolumeCounter_Y_First = 0;
#if defined(DEBUG_15)
    VolumeCounter_Y_Last  = 0;
#else
    VolumeCounter_Y_Last  = Volume_NumY - 1;
#endif
    VolumeCounter_Y_IncrementDirection =  1;

        Volume_Mu_XZSlab_Current_ptr =     Volume_Mu_Active_ptr;
    Volume_Energy_XZSlab_Current_ptr = Volume_Energy_Active_ptr;
    }
  else
    {
    VolumeCounter_Y_First = Volume_NumY - 1;
#if defined(DEBUG_15)
    VolumeCounter_Y_Last  = Volume_NumY - 1;
#else
    VolumeCounter_Y_Last  = 0;
#endif
    VolumeCounter_Y_IncrementDirection = -1;

        Volume_Mu_XZSlab_Current_ptr =     Volume_Mu_Active_ptr + Volume_NumVoxels_XZSlab*(Volume_NumY - 1);
    Volume_Energy_XZSlab_Current_ptr = Volume_Energy_Active_ptr + Volume_NumVoxels_XZSlab*(Volume_NumY - 1);
    }

  // Done initializing direction to progress in the Y direction.
  // *******************************************************************************************************************
  // *******************************************************************************************************************
  // For each row in the image volume (i.e. XZ slab at constant Y), project the X and Z coordinates of the image volume
  // voxel boundaries onto the XZ plane (at Y = 0), and project the incident energy view into that voxel slab.
  
#if defined(DEBUG_00)
  Report("For each row in the image volume (i.e. XZ slab at constant Y), project the X and Z coordinates of the image volume\n");
  Report("voxel boundaries onto the XZ plane (at Y = 0), and backproject the current 2D view into that voxel slab.\n");
#endif

  VolumeCounter_Y_BeyondLast  = VolumeCounter_Y_Last + VolumeCounter_Y_IncrementDirection;
  for (VolumeCounter_Y  = VolumeCounter_Y_First;
       VolumeCounter_Y != VolumeCounter_Y_BeyondLast;
       VolumeCounter_Y += VolumeCounter_Y_IncrementDirection)
    {
    // The position of the center of the current XZ slab relative to the volume center:
                     VolumeXZSlab_Y_Current = Volume_VoxelSize_Y*(Volume_NumY/2 - VolumeCounter_Y - 1/2.);
    VolumeXZSlabMagnificationFactor_Current = SourceCoord_Y/(SourceCoord_Y - VolumeXZSlab_Y_Current);
    
    // Initialize image parameters.
    VolumeVoxelBoundaryCoordProjected_X_Step    = Volume_VoxelSize_X * VolumeXZSlabMagnificationFactor_Current;
    VolumeVoxelBoundaryCoordProjected_Z_Step    = Volume_VoxelSize_Z * VolumeXZSlabMagnificationFactor_Current;
    
    VolumeVoxelBoundaryCoordProjected_X_Initial = SourceCoord_X - (Volume_NumX/2. + SourceCoord_X)*VolumeVoxelBoundaryCoordProjected_X_Step;
    
    VolumeVoxelBoundaryCoordProjected_Z_Initial = SourceCoord_Z * (VolumeXZSlabMagnificationFactor_Current - 1)
                                                - VolumeVoxelBoundaryCoordProjected_Z_Step * Volume_NumZ/2.;
    
    // Detector_AbsorbedEnergyViewData needs to be cleared for each slab. Note that the un-reversed (not _Current) pointer is correct.
    memset(Detector_AbsorbedEnergyViewData_ptr, 0, Detector_NumPixels*sizeof(float));

#if defined(DEBUG_16)
    sprintf(OutputString, "Backprojecting a volume XZ slab at VolumeCounter_Y = % 3i, from SourceCoord_X = % 8.2f, SourceCoord_Y = % 8.2f, SourceCoord_Z = % 8.2f\n", 
                                                                VolumeCounter_Y,             SourceCoord_X,          SourceCoord_Y,          SourceCoord_Z);
    Report(OutputString);
#endif
#if defined(DEBUG_17)
    sprintf(OutputString, "Detector_NumCols = %3i, Detector_NumRows = %3i, Detector_IncrementDirection = % 3i\n",
                           Detector_NumCols,       Detector_NumRows,       Detector_IncrementDirection);
    Report(OutputString);
#endif
#if defined(DEBUG_18)
    sprintf(OutputString, "VolumeCounter_Y_First = %3i, VolumeCounter_Y_Last = %3i, VolumeCounter_Y_IncrementDirection = % 3i\n",
                           VolumeCounter_Y_First,       VolumeCounter_Y_Last,       VolumeCounter_Y_IncrementDirection);
    Report(OutputString);
    sprintf(OutputString, "Volume_VoxelSize_X = % 8.2f, Volume_VoxelSize_Y = % 8.2f, Volume_VoxelSize_Z = % 8.2f\n",
                           Volume_VoxelSize_X,          Volume_VoxelSize_Y,          Volume_VoxelSize_Z);
    Report(OutputString);
    sprintf(OutputString, "VolumeXZSlab_Y_Current = % 8.2f, VolumeXZSlabMagnificationFactor_Current = % 8.2f\n",
                           VolumeXZSlab_Y_Current,          VolumeXZSlabMagnificationFactor_Current);
    Report(OutputString);
#endif
#if defined(DEBUG_19)
    sprintf(OutputString, "DetectorPixelBoundaryCoordRotated_X            = % 8.2f, DetectorPixelBoundaryCoordRotated_X (center)    = % 8.2f, DetectorPixelBoundaryCoordRotated_X (last)    = % 8.2f\n",
                          *(DetectorPixelBoundaryCoordRotated_X_Current_ptr                                                        ),
                          *(DetectorPixelBoundaryCoordRotated_X_Current_ptr    + Detector_IncrementDirection*Detector_NumCols/2    ),
                          *(DetectorPixelBoundaryCoordRotated_X_Current_ptr    + Detector_IncrementDirection*Detector_NumCols      ));
    Report(OutputString);
    sprintf(OutputString, "DetectorPixelBoundaryCoordRotated_Y            = % 8.2f, DetectorPixelBoundaryCoordRotated_Y (center)    = % 8.2f, DetectorPixelBoundaryCoordRotated_Y (last)    = % 8.2f\n",
                          *(DetectorPixelBoundaryCoordRotated_Y_Current_ptr                                                        ),
                          *(DetectorPixelBoundaryCoordRotated_Y_Current_ptr    + Detector_IncrementDirection*Detector_NumCols/2    ),
                          *(DetectorPixelBoundaryCoordRotated_Y_Current_ptr    + Detector_IncrementDirection*Detector_NumCols      ));
    Report(OutputString);
    sprintf(OutputString, "DetectorPixelBoundaryCoordsShifted_Z           = % 8.2f, DetectorPixelBoundaryCoordsShifted_Z (center)   = % 8.2f, DetectorPixelBoundaryCoordsShifted_Z (last)   = % 8.2f\n",
                          *(DetectorPixelBoundaryCoordsShifted_Z_ptr                                                      ),
                          *(DetectorPixelBoundaryCoordsShifted_Z_ptr           +                    Detector_NumRows/2    ),
                          *(DetectorPixelBoundaryCoordsShifted_Z_ptr           +                    Detector_NumRows      ));
    Report(OutputString);
    sprintf(OutputString, "DetectorMagnificationFactor_Current            = % 8.2f, DetectorMagnificationFactor (center-ish)        = % 8.2f, DetectorMagnificationFactor (last)            = % 8.2f, DetectorMagnificationFactor (sentinal)            = % 8.2f\n",
                          *(DetectorMagnificationFactor_Current_ptr                                                                ),
                          *(DetectorMagnificationFactor_Current_ptr            + Detector_IncrementDirection*Detector_NumCols/2    ),
                          *(DetectorMagnificationFactor_Current_ptr            + Detector_IncrementDirection*Detector_NumCols      ),
                          *(DetectorMagnificationFactor_Current_ptr            + Detector_IncrementDirection*Detector_NumCols   + 1));
    Report(OutputString);
    sprintf(OutputString, "DetectorPixelBoundaryCoordsProjected_X_Initial = % 8.2f, DetectorPixelBoundaryCoordsProjected_X (center) = % 8.2f, DetectorPixelBoundaryCoordsProjected_X (last) = % 8.2f, DetectorPixelBoundaryCoordsProjected_X (sentinal) = % 8.2f\n",
                          *(DetectorPixelBoundaryCoordsProjected_X_Initial_ptr                                                     ),
                          *(DetectorPixelBoundaryCoordsProjected_X_Initial_ptr + Detector_IncrementDirection*Detector_NumCols/2 + 1),
                          *(DetectorPixelBoundaryCoordsProjected_X_Initial_ptr + Detector_IncrementDirection*Detector_NumCols      ),
                          *(DetectorPixelBoundaryCoordsProjected_X_Initial_ptr + Detector_IncrementDirection*Detector_NumCols   + 1));
    Report(OutputString);
#endif
#if defined(DEBUG_20)
    sprintf(OutputString, "VolumeVoxelBoundaryCoordProjected_X_Step = % 8.2f, Volume_NumX = % 3i\n",
                          VolumeVoxelBoundaryCoordProjected_X_Step, Volume_NumX);
    Report(OutputString);
    sprintf(OutputString, "VolumeVoxelBoundaryCoordProjected_Z_Step = % 8.2f, Volume_NumZ = % 3i\n",
                          VolumeVoxelBoundaryCoordProjected_Z_Step, Volume_NumZ);
    Report(OutputString);
    sprintf(OutputString, "VolumeVoxelBoundaryCoordProjected_X_Initial    = % 8.2f, VolumeVoxelBoundaryCoordProjected_X (center)    = % 8.2f, VolumeVoxelBoundaryCoordProjected_X (last)    = % 8.2f\n",
                          VolumeVoxelBoundaryCoordProjected_X_Initial,
                          VolumeVoxelBoundaryCoordProjected_X_Initial + VolumeVoxelBoundaryCoordProjected_X_Step*Volume_NumX/2,
                          VolumeVoxelBoundaryCoordProjected_X_Initial + VolumeVoxelBoundaryCoordProjected_X_Step*Volume_NumX  );
    Report(OutputString);
    sprintf(OutputString, "VolumeVoxelBoundaryCoordProjected_Z_Initial    = % 8.2f, VolumeVoxelBoundaryCoordProjected_Z (center)    = % 8.2f, VolumeVoxelBoundaryCoordProjected_Z (last)    = % 8.2f\n",
                          VolumeVoxelBoundaryCoordProjected_Z_Initial,
                          VolumeVoxelBoundaryCoordProjected_Z_Initial + VolumeVoxelBoundaryCoordProjected_Z_Step*Volume_NumZ/2,
                          VolumeVoxelBoundaryCoordProjected_Z_Initial + VolumeVoxelBoundaryCoordProjected_Z_Step*Volume_NumZ  );
    Report(OutputString);
#endif
#if defined(DEBUG_21)
    sprintf(OutputString, "Detector_IncidentEnergyViewData (center-ish) = % 8.3g, Volume_Mu_XZSlab (center-ish) = % 8.3g, Volume_Energy_XZSlab (center-ish) (before projection of this XZ slab) = % 8.3g\n",
                          *(Detector_IncidentEnergyViewData_ptr + Detector_IncrementDirection*Detector_NumPixels/2 + Detector_IncrementDirection*Detector_NumRows/2),
    //                      beginning of projection             + halfway through the projection (to the beginning of a detector column)
    //                                                                                                             + halfway across along that column
                          *(Volume_Mu_XZSlab_Current_ptr        + Volume_NumVoxels_XZSlab/2 + Volume_NumZ/2),
    //                      beginning of XZ slab                + halfway through the slab  (to the beginning of the 1st Z-stack of voxels in the 2nd half of the slab)
    //                                                                                      + halfway along that Z-stack
                          *(Volume_Energy_XZSlab_Current_ptr    + Volume_NumVoxels_XZSlab/2 + Volume_NumZ/2));
    //                      same arithmetic as Volume_Mu_XZSlab_Current
    Report(OutputString);
#endif

#if defined(DEBUG_51)
    sprintf(OutputString, "\nBefore projecting a volume XZ slab at VolumeCounter_Y = % 3i\n", VolumeCounter_Y);
    Report(OutputString);
    PrintArray("\nIncident energy projection:\n",   Detector_IncidentEnergyViewData_Current_ptr, 1e-4, 0, Detector_NumCols, Detector_NumRows, Detector_IncrementDirection);
    PrintArray("\nAbsorbed energy projection:\n",   Detector_AbsorbedEnergyViewData_Current_ptr, 0.0,  0, Detector_NumCols, Detector_NumRows, Detector_IncrementDirection);
    PrintArray("\nMu Volume XZ slab:\n",            Volume_Mu_XZSlab_Current_ptr,                0.0,  0, Volume_NumX,      Volume_NumZ,      1);
    PrintArray("\nEnergy Volume XZ slab:\n",        Volume_Energy_XZSlab_Current_ptr,            0.0,  0, Volume_NumX,      Volume_NumZ,      1);
#endif

    DD3EnergyRow(VolumeVoxelBoundaryCoordProjected_X_Initial, VolumeVoxelBoundaryCoordProjected_X_Step, Volume_NumX,
                 VolumeVoxelBoundaryCoordProjected_Z_Initial, VolumeVoxelBoundaryCoordProjected_Z_Step, Volume_NumZ,
                 Volume_Mu_XZSlab_Current_ptr, Volume_Energy_XZSlab_Current_ptr,
                 DetectorPixelBoundaryCoordsProjected_X_Initial_ptr, Detector_IncrementDirection,
                 DetectorPixelBoundaryCoordsShifted_Z_ptr, DetectorMagnificationFactor_Current_ptr, SourceCoord_Z,
                 Detector_IncidentEnergyViewData_Current_ptr, Detector_AbsorbedEnergyViewData_Current_ptr, Detector_PathLength_Current_ptr,
                 Detector_NumCols, Detector_NumRows);
    
#if defined(DEBUG_52)
    sprintf(OutputString, "\n After projecting a volume XZ slab at VolumeCounter_Y = % 3i\n", VolumeCounter_Y);
    Report(OutputString);
    PrintArray("\nIncident energy projection:\n",   Detector_IncidentEnergyViewData_Current_ptr, 1e-4, 0, Detector_NumCols, Detector_NumRows, Detector_IncrementDirection);
    PrintArray("\nAbsorbed energy projection:\n",   Detector_AbsorbedEnergyViewData_Current_ptr, 0.0,  0, Detector_NumCols, Detector_NumRows, Detector_IncrementDirection);
    PrintArray("\nMu Volume XZ slab:\n",            Volume_Mu_XZSlab_Current_ptr,                0.0,  0, Volume_NumX,      Volume_NumZ,      1);
    PrintArray("\nEnergy Volume XZ slab:\n",        Volume_Energy_XZSlab_Current_ptr,            0.0,  0, Volume_NumX,      Volume_NumZ,      1);
#endif

#if defined(DEBUG_23)
    sprintf(OutputString, "Detector_IncidentEnergyViewData (center-ish) = % 8.3g, Volume_Mu_XZSlab (center-ish) = % 8.3g, Volume_Energy_XZSlab (center-ish)  (after projection of this XZ slab) = % 8.3g\n",
                          *(Detector_IncidentEnergyViewData_ptr + Detector_IncrementDirection*Detector_NumPixels/2 + Detector_IncrementDirection*Detector_NumRows/2),
                          *(Volume_Mu_XZSlab_Current_ptr        + Volume_NumVoxels_XZSlab/2 + Volume_NumZ/2),
                          *(Volume_Energy_XZSlab_Current_ptr    + Volume_NumVoxels_XZSlab/2 + Volume_NumZ/2));
    Report(OutputString);
#endif

#if defined(DEBUG_24)
    if ((VolumeCounter_Y == 0)
     || (VolumeCounter_Y == (1*Volume_NumY/10))
     || (VolumeCounter_Y == (2*Volume_NumY/10))
     || (VolumeCounter_Y == (3*Volume_NumY/10))
     || (VolumeCounter_Y == (4*Volume_NumY/10))
     || (VolumeCounter_Y == (5*Volume_NumY/10))
     || (VolumeCounter_Y == (6*Volume_NumY/10))
     || (VolumeCounter_Y == (7*Volume_NumY/10))
     || (VolumeCounter_Y == (8*Volume_NumY/10))
     || (VolumeCounter_Y == (9*Volume_NumY/10))
     || (VolumeCounter_Y == VolumeCounter_Y_Last))
      {
      sprintf(OutputString, "\n After projecting a volume XZ slab at VolumeCounter_Y = % 3i\n", VolumeCounter_Y);
      Report(OutputString);
      
      int Counter_X, Counter_Y, Counter_Z;
      
      for (Counter_Z=Volume_NumZ/2 ; Counter_Z <= (Volume_NumZ-1) ; Counter_Z+=Volume_NumZ)
        {
        sprintf(OutputString, "Slice % 3i\n", Counter_Z);
        Report(OutputString);
        for (Counter_Y=0 ; Counter_Y <= (Volume_NumY-1) ; Counter_Y+=Volume_NumY/10)
          {
          for (Counter_X=0 ; Counter_X <= (Volume_NumX-1) ; Counter_X+=Volume_NumX/10)
            {
            sprintf(OutputString, "%10.5f ", *(Volume_Energy_Active_ptr + Volume_NumVoxels_XZSlab*Counter_Y + Volume_NumZ*Counter_X + Counter_Z));
            Report(OutputString);
            }
          sprintf(OutputString, "\n");
          Report(OutputString);
          }
        }
      }
#endif

    // Done projecting XZ slab.
    // *******************************************************************************************************************
    // *******************************************************************************************************************
    // Increment to the next XZ slab in the volume.

#if defined(DEBUG_00)
    Report("Increment to the next XZ slab in the volume.\n");
#endif
    Volume_Mu_XZSlab_Current_ptr     += VolumeCounter_Y_IncrementDirection * Volume_NumVoxels_XZSlab;
    Volume_Energy_XZSlab_Current_ptr += VolumeCounter_Y_IncrementDirection * Volume_NumVoxels_XZSlab;
    
    } // end volume Y (i.e. XZ-slab) loop
  // Done projecting the view onto each XZ slab in the volume.
  // *******************************************************************************************************************
  
  }

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************

// This function performs a distance-driven ptojection of a 3D sinogram of incident energy into a 3D volume of Mu values,
// to produce a 3D volume of deposited energies.
// NOTE: View-by-view depositied energy is SUBTRACTED from the input sinogram containing incident energy. Therefore, on
// exit, the input sinogram contains the energy that was NOT absorbed, therefore this should be DETECTED energy.

extern "C" {
EXPORT
void DD3Energy(float  SourceCoord_X,                               // Input:  scalar
	             float  SourceCoord_Y,                               // Input:  scalar
	             float  SourceCoord_Z,                               // Input:  scalar
	             int    Detector_NumCols,                            // Input:  scalar
	             int    Detector_NumRows,                            // Input:  scalar
	             float *DetectorPixelCenterCoords_X_ptr,             // Input:  [Detector_NumCols]
	             float *DetectorPixelCenterCoords_Y_ptr,             // Input:  [Detector_NumCols]
	             float *DetectorPixelCenterCoords_Z_ptr,             // Input:  [Detector_NumRows]
	             float  DetectorPixelAspectRatio_Z_X,                // Input:  scalar
	             float  VolumeOffset_X,                              // Input:  scalar
	             float  VolumeOffset_Y,                              // Input:  scalar
	             float  VolumeOffset_Z,                              // Input:  scalar
  	           float *ViewAngles_ptr,                              // Input:  [NumViews]
	             float *Shifts_Z_EachView_ptr,                       // Input:  [NumViews]
	             int    NumViews,                                    // Input:  scalar
	             float *Detector_IncidentEnergySinogramData_ptr,     // Input:  [NumViews][Detector_NumCols][Detector_NumRows]
	             int    Volume_NumX,                                 // Input:  scalar
	             int    Volume_NumY,                                 // Input:  scalar
	             int    Volume_NumZ,                                 // Input:  scalar
	             float  Volume_VoxelSize_X,                          // Input:  scalar
	             float  Volume_VoxelSize_Y,                          // Input:  scalar
	             float  Volume_VoxelSize_Z,                          // Input:  scalar
	             float *Volume_Mu_Untransposed_ptr,                  // Input:  [Volume_NumY][Volume_NumX][Volume_NumZ]
	             float *Volume_Energy_Untransposed_ptr)              // Output: [Volume_NumY][Volume_NumX][Volume_NumZ]
  {
  int    Detector_NumPixels;
  int    Detector_NumPixelBoundaries_XY;
  int    Detector_NumPixelBoundaries_Z;
  int    Counter;
  int    ViewCounter;
  int    ViewCounter_Last;
  int    ProjectionVertical_flag;
  float  sinViewAngle;
  float  cosViewAngle;
  float  SourceCoordRotated_X;
  float  SourceCoordRotated_Y;
  float  SourceCoordShifted_Z;
  float *DetectorPixelBoundaryCoords_X_ptr;
  float *DetectorPixelBoundaryCoords_Y_ptr;
  float *DetectorPixelBoundaryCoords_Z_ptr;
  float *DetectorPixelBoundaryCoord_X_Current_ptr;
  float *DetectorPixelBoundaryCoord_Y_Current_ptr;
  float *DetectorPixelBoundaryCoord_Z_Current_ptr;
  float *DetectorPixelBoundaryCoordsRotated_X_ptr;
  float *DetectorPixelBoundaryCoordsRotated_Y_ptr;
  float *DetectorPixelBoundaryCoordsShifted_Z_ptr;
  float *DetectorPixelBoundaryCoordRotated_X_Current_ptr;
  float *DetectorPixelBoundaryCoordRotated_Y_Current_ptr;
  float *DetectorPixelBoundaryCoordShifted_Z_Current_ptr;
  float *DetectorPixelBoundaryCoordsProjected_X_ptr;
  float *DetectorMagnificationFactors_ptr;
  float *ViewAngle_Current_ptr;
  float *Detector_PathLengths_ptr;
  float *Detector_IncidentEnergyViewData_Current_ptr;
  float *Detector_AbsorbedEnergyViewData_Current_ptr;
  float *Volume_Mu_Transposed_ptr;
  float *Volume_Energy_Transposed_ptr;

#if defined(DEBUG_00)
  sprintf(OutputString, "Running DD3Energy\n");
  Report(OutputString);
#endif

#if defined(DEBUG_10)
  int   ViewCounter_First;

  ViewCounter_First = 0;
  ViewCounter_Last  = 0;
#else
  ViewCounter_Last = NumViews - 1;
#endif

  Detector_NumPixels             = Detector_NumCols * Detector_NumRows;
  Detector_NumPixelBoundaries_XY = Detector_NumCols+1;
  Detector_NumPixelBoundaries_Z  = Detector_NumRows+1;

  // Allocate memory for detector boundaries
  DetectorPixelBoundaryCoords_X_ptr          = (float*)malloc((Detector_NumPixelBoundaries_XY    ) * sizeof(float));
  DetectorPixelBoundaryCoords_Y_ptr          = (float*)malloc((Detector_NumPixelBoundaries_XY    ) * sizeof(float));
  DetectorPixelBoundaryCoords_Z_ptr          = (float*)malloc((Detector_NumPixelBoundaries_Z     ) * sizeof(float));
  DetectorPixelBoundaryCoordsRotated_X_ptr   = (float*)malloc((Detector_NumPixelBoundaries_XY    ) * sizeof(float));
  DetectorPixelBoundaryCoordsRotated_Y_ptr   = (float*)malloc((Detector_NumPixelBoundaries_XY    ) * sizeof(float));
  DetectorPixelBoundaryCoordsShifted_Z_ptr   = (float*)malloc((Detector_NumPixelBoundaries_Z  + 2) * sizeof(float)); // 2 sentinels
  
  // Allocate memory for transposed volumes.
      Volume_Mu_Transposed_ptr = (float*)calloc(Volume_NumX*Volume_NumY*Volume_NumZ,sizeof(float));
  Volume_Energy_Transposed_ptr = (float*)calloc(Volume_NumX*Volume_NumY*Volume_NumZ,sizeof(float));

  // The "vertical" projections need Volume_Mu_Untransposed (the input), whereas
  // the "horizontal" projections need Volume_Mu_Transposed.
  // So we need to make a transposed copy of Volume_Mu_Untransposed.
  DD3AddTranspose(Volume_NumY, Volume_NumX, Volume_NumZ, Volume_Mu_Untransposed_ptr, Volume_Mu_Transposed_ptr);

  // The next four arrays are passed to a called function for use there.
  // The could be allocated in that function, but they are allocated here for efficiency.
  
             DetectorMagnificationFactors_ptr = (float*)malloc((Detector_NumCols               + 2) * sizeof(float)); // 2 sentinels
   DetectorPixelBoundaryCoordsProjected_X_ptr = (float*)malloc((Detector_NumPixelBoundaries_XY + 2) * sizeof(float)); // 2 sentinels
                     Detector_PathLengths_ptr = (float*)calloc((Detector_NumPixels                ),  sizeof(float));
  Detector_AbsorbedEnergyViewData_Current_ptr = (float*)calloc((Detector_NumPixels                ),  sizeof(float));

  // Calculate detector boundaries
  DD3Boundaries(Detector_NumPixelBoundaries_XY, DetectorPixelCenterCoords_X_ptr, DetectorPixelBoundaryCoords_X_ptr);
  DD3Boundaries(Detector_NumPixelBoundaries_XY, DetectorPixelCenterCoords_Y_ptr, DetectorPixelBoundaryCoords_Y_ptr);
  DD3Boundaries(Detector_NumPixelBoundaries_Z,  DetectorPixelCenterCoords_Z_ptr, DetectorPixelBoundaryCoords_Z_ptr);

  // Translate the Z position of each boundary by the Z position of the source.
  for (Counter=0 ; Counter <= (Detector_NumPixelBoundaries_Z-1) ; Counter++)
    *(DetectorPixelBoundaryCoords_Z_ptr + Counter) -= SourceCoord_Z;
  
  // Translate the Z position of the source by the specified Z offset.
  SourceCoord_Z -= VolumeOffset_Z;

  // *******************************************************************************************************************
  // Loop through all views.
  
#if defined(DEBUG_00)
  Report("Loop through all views.\n");
#endif

  Detector_IncidentEnergyViewData_Current_ptr = Detector_IncidentEnergySinogramData_ptr;
                        ViewAngle_Current_ptr =                          ViewAngles_ptr;

  for (ViewCounter = 0 ; ViewCounter <= ViewCounter_Last ; ViewCounter++)
    {
#if defined(DEBUG_10)
    if ((ViewCounter == 0) && (ViewCounter_First > 0))
      {
      Detector_IncidentEnergyViewData_Current_ptr += (ViewCounter_First * Detector_NumPixels);
                            ViewAngle_Current_ptr +=  ViewCounter_First;
                            Shifts_Z_EachView_ptr +=  ViewCounter_First;
                            ViewCounter            =  ViewCounter_First;
      }

    if (ViewCounter == ViewCounter_Last)
      PrintDebug = 1;
    else
      PrintDebug = 0;
#endif

    sinViewAngle = (float)sin(*ViewAngle_Current_ptr);
    cosViewAngle = (float)cos(*ViewAngle_Current_ptr);
    
    DetectorPixelBoundaryCoord_X_Current_ptr        = DetectorPixelBoundaryCoords_X_ptr;
    DetectorPixelBoundaryCoord_Y_Current_ptr        = DetectorPixelBoundaryCoords_Y_ptr;
    DetectorPixelBoundaryCoordRotated_X_Current_ptr = DetectorPixelBoundaryCoordsRotated_X_ptr;
    DetectorPixelBoundaryCoordRotated_Y_Current_ptr = DetectorPixelBoundaryCoordsRotated_Y_ptr;
    
    for (Counter=0 ; Counter<=(Detector_NumPixelBoundaries_XY-1) ; Counter++)
      {
      // Rotate detector pixel X coordinate.
      *DetectorPixelBoundaryCoordRotated_X_Current_ptr++ = *DetectorPixelBoundaryCoord_X_Current_ptr   * cosViewAngle
                                                         - *DetectorPixelBoundaryCoord_Y_Current_ptr   * sinViewAngle
                                                         - VolumeOffset_X;
      // Rotate detector pixel Y coordinates.
      *DetectorPixelBoundaryCoordRotated_Y_Current_ptr++ = *DetectorPixelBoundaryCoord_Y_Current_ptr++ * cosViewAngle
                                                         + *DetectorPixelBoundaryCoord_X_Current_ptr++ * sinViewAngle
                                                         - VolumeOffset_Y;
      }
    
    SourceCoordRotated_X    = SourceCoord_X * cosViewAngle - SourceCoord_Y * sinViewAngle;
    SourceCoordRotated_Y    = SourceCoord_Y * cosViewAngle + SourceCoord_X * sinViewAngle;
    ProjectionVertical_flag = (fabs(SourceCoordRotated_Y) >= fabs(SourceCoordRotated_X));
    SourceCoordRotated_X   -= VolumeOffset_X;
    SourceCoordRotated_Y   -= VolumeOffset_Y;

    // Shift z coordinates
    DetectorPixelBoundaryCoord_Z_Current_ptr        = DetectorPixelBoundaryCoords_Z_ptr;
    DetectorPixelBoundaryCoordShifted_Z_Current_ptr = DetectorPixelBoundaryCoordsShifted_Z_ptr;

    for (Counter=0 ; Counter<=(Detector_NumPixelBoundaries_Z-1) ; Counter++)
      *DetectorPixelBoundaryCoordShifted_Z_Current_ptr++ = *DetectorPixelBoundaryCoord_Z_Current_ptr++; // + *Shifts_Z_EachView_ptr;

        *(DetectorPixelBoundaryCoordsShifted_Z_ptr + Detector_NumPixelBoundaries_Z    )
      = *(DetectorPixelBoundaryCoordsShifted_Z_ptr + Detector_NumPixelBoundaries_Z - 1); 
//     *(DetectorPixelBoundaryCoordsShifted_Z_ptr + Detector_NumPixelBoundaries_Z    ) = 1.0e12;  // sentinel in shifted array 
//     *(DetectorPixelBoundaryCoordsShifted_Z_ptr + Detector_NumPixelBoundaries_Z + 1) = 1.5e12;  // sentinel in shifted array 
    SourceCoordShifted_Z = SourceCoord_Z + *Shifts_Z_EachView_ptr;

#if defined(DEBUG_05) || defined(DEBUG_16) || defined(DEBUG_17) || defined(DEBUG_18) || defined(DEBUG_19) || defined(DEBUG_20) || defined(DEBUG_21) || defined(DEBUG_22) || defined(DEBUG_23) || defined(DEBUG_24) 
    sprintf(OutputString, "\nBackprojecting a view: ViewCounter = % 3i, ViewAngle_Current = % 8.2f (% 3.1f degrees), ProjectionVertical_flag = % 3i\n",
                          ViewCounter, *ViewAngle_Current_ptr, *ViewAngle_Current_ptr*180./3.141592, ProjectionVertical_flag);
    Report(OutputString);
#endif

    DD3EnergyView(SourceCoordRotated_X, SourceCoordRotated_Y, SourceCoordShifted_Z,
                  Detector_NumCols, Detector_NumRows,
                  ProjectionVertical_flag,
                  DetectorPixelBoundaryCoordsRotated_X_ptr,
                  DetectorPixelBoundaryCoordsRotated_Y_ptr,
                  DetectorPixelBoundaryCoordsProjected_X_ptr,
                  DetectorPixelBoundaryCoordsShifted_Z_ptr,
                  DetectorMagnificationFactors_ptr, DetectorPixelAspectRatio_Z_X,
                  Detector_IncidentEnergyViewData_Current_ptr, Detector_AbsorbedEnergyViewData_Current_ptr, Detector_PathLengths_ptr,
                  Volume_NumX, Volume_NumY, Volume_NumZ,
                  Volume_VoxelSize_X, Volume_VoxelSize_Y, Volume_VoxelSize_Z,
                  Volume_Mu_Untransposed_ptr, Volume_Mu_Transposed_ptr,
                  Volume_Energy_Untransposed_ptr, Volume_Energy_Transposed_ptr);
    
    // Move to the next view.
    Detector_IncidentEnergyViewData_Current_ptr += Detector_NumPixels;
    ++Shifts_Z_EachView_ptr;
    ++ViewAngle_Current_ptr;
    
    }            
  // end of view loop

  // Results from the "vertical" backprojections got written into Volume_Energy_Untransposed, whereas
  // results from the "horizontal" backprojections got written into Volume_Energy_Transposed.
  // So we need to transpose the transposed volume, and add that to the un-transposed volume.
  DD3AddTranspose(Volume_NumY, Volume_NumX, Volume_NumZ, Volume_Energy_Transposed_ptr, Volume_Energy_Untransposed_ptr);

  // Clean up memory
  free(DetectorPixelBoundaryCoords_X_ptr);
  free(DetectorPixelBoundaryCoords_Y_ptr);
  free(DetectorPixelBoundaryCoords_Z_ptr);
  free(DetectorPixelBoundaryCoordsRotated_X_ptr);
  free(DetectorPixelBoundaryCoordsRotated_Y_ptr);
  free(DetectorPixelBoundaryCoordsShifted_Z_ptr);
  free(Volume_Mu_Transposed_ptr);
  free(Volume_Energy_Transposed_ptr);
  free(DetectorMagnificationFactors_ptr);
  free(DetectorPixelBoundaryCoordsProjected_X_ptr);
  free(Detector_PathLengths_ptr);
  free(Detector_AbsorbedEnergyViewData_Current_ptr);
  
  } // end of function
} // end of extern

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
