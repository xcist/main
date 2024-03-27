// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef __DD3_roi_notrans_mm__
#define __DD3_roi_notrans_mm__

#define MIN_(a,b) (((a) < (b)) ? (a) : (b))
#define MAX_(a,b) (((a) > (b)) ? (a) : (b))

extern bool useUInt16;

/*
* DD3 transpose
*/
void DD3Transpose(int nrcols,
				  int nrrows,
				  int nrplanes,
				  void *pOrig,
				  float *pTrans);

/*
* DD3 boundaries
*/
void DD3Boundaries(int nrBoundaries,
				   float *pCenters,
				   float *pBoundaries);

typedef unsigned char byte;


//roi: only apply on a part of the image
//notrans: no pTans for 3D volume projector and backprojector
//mm: all the metric is in mm


void DD3ProjRow_roi_notrans_mm(float imgX,
								float imgXstep,
								int nrcols,
								float imgZstart,
								float imgZstep,
								int nrplanes,
								void *pImg,
								float *detX,
								int increment,
								float *detZ,
								float *scales,
								float z0,
								float *view,
								int nrdetrows,
								int nrdetcols,
								int trans_rowstep,
								byte* xy_mask);

void DD3ProjView_roi_notrans_mm(float x0,
								 float y0,
								 float z0,
								 int nrdetcols,
								 int nrdetrows,
								 int vertical,
								 float *xdi,
								 float *ydi,
								 float *detX, // empty distance array (nrxdist + 2)
								 float *detZ,
								 float *scales, // empty array (nrdetcols + 2)
								 float *sinogram,
								 float *view, // empty view with 1 pixel margin all around
								 int nrcols,
								 int nrrows,       // images
								 int nrplanes,     //     do NOT
								 void *pOrig,     //        contain
								 float vox_xy_size,   // voxel size
								 float vox_z_size,
								 byte* xy_mask,
								 byte* xy_mask_trans);    //field to represent vox sizes

void DD3BackRow_roi_notrans_mm(float imgX,
								float imgXstep,
								int nrcols,
								float imgZstart,
								float imgZstep,
								int nrplanes,
								void *pImg,
								float *detX,
								int increment,
								float *detZ,
								float *scales,
								float z0,
								float *view,
								int nrdetrows,
								int nrdetcols,
								int trans_rowstep,
								byte* xy_mask);


void DD3BackView_roi_notrans_mm(float x0,
								 float y0,
								 float z0,
								 int nrdetcols,
								 int nrdetrows,
								 int vertical,
								 float *xdi,
								 float *ydi,
								 float *detX, // empty distance array (nrxdist + 2)
								 float *detZ,
								 float *scales, // empty array (nrdetcols + 2)
								 float *sinogram,
								 float *view, // empty view with 1 pixel margin all around
								 int nrcols,
								 int nrrows,       // images
								 int nrplanes,     //     do NOT
								 void *pOrig,     //        contain
								 float vox_xy_size,   // voxel size
								 float vox_z_size,
								 byte* xy_mask,
								 byte* xy_mask_trans);

typedef struct
{
	float x0;
	float y0;
	float z0;
	int nrdetcols;
	int nrdetrows;
	float *xds;
	float *yds;
	float *zds;
	float *viewangles;
    float *zshifts;
	int nrviews;
}DD3ScanGeo; //Read only

typedef struct  
{
	float imgXoffset;
	float imgYoffset;
	float imgZoffset;
	int nrcols;
	int nrrows;
	int nrplanes;
	float vox_xy_size;
	float vox_z_size;
	byte *xy_mask;
}DD3ImgGeo; //Read only

typedef struct
{
	float *xdi; /* ns+1 */
	float *ydi; /* ns+1 */
	float *xdiRot; /* ns+1 */
	float *ydiRot; /* ns+1 */
	float *detZ; /* nt+1+1 */
	float *view; /* (ns+2)*(nt+2) */
	float *detX; /* ns+3 */
	float *scales; /* ns+2 */
} DD3GeoBuff; //Buffer to modify


void DD3SetBuff(DD3ScanGeo &ScanGeo, DD3ImgGeo &ImgGeo, DD3GeoBuff &GeoBuff);
void DD3Proj_struct_kernel(DD3ScanGeo ScanGeo, DD3ImgGeo ImgGeo, DD3GeoBuff GeoBuff,
						   float* sinogram,float* img);

void DD3Back_struct_kernel(DD3ScanGeo ScanGeo, DD3ImgGeo ImgGeo, DD3GeoBuff GeoBuff,
						   float* sinogram,float* pOrig);

void DD3FreeBuff(DD3GeoBuff GeoBuf);


/*
* DD3 projector and backprojector
*/
extern "C"{

	void DD3Proj_roi_notrans_mm(float x0,
		float y0,
		float z0,
		int nrdetcols,
		int nrdetrows,
		float *xds,
		float *yds,
		float *zds,
		float imgXoffset,
		float imgYoffset,
		float imgZoffset,
		float *viewangles,
		float *zshifts,
		int nrviews,
		float *sinogram,
		int nrcols,         // image
		int nrrows,         //    does NOT
		int nrplanes,       //        contain a dummy 1 pixel frame
		void *pOrig,
		float vox_xy_size, //added fields
		float vox_z_size,//added fields 
		byte* xy_mask); //added fields for xy_mask
		//float* xy_wob);//HD support for focal spot wob

	void DD3Back_roi_notrans_mm(float x0,
		float y0,
		float z0,
		int nrdetcols,
		int nrdetrows,
		float *xds,
		float *yds,
		float *zds,
		float imgXoffset,
		float imgYoffset,
		float imgZoffset,
		float *viewangles,
		float *zshifts,
		int nrviews,
		float *sinogram,
		int nrcols,         // image
		int nrrows,         //    does NOT
		int nrplanes,       //        contain a dummy 1 pixel frame
		void *pOrig,
		float vox_xy_size, //added fields
		float vox_z_size,//added fields 
		byte* xy_mask); //added fields for xy_mask
		//float *xy_wob); //HD support for focal spot wob


	//////////////////////////////////////////////////////////////////////////
	//structure parameters callback functions
	void DD3Proj_struct(float x0,
		float y0,
		float z0,
		int nrdetcols,
		int nrdetrows,
		float *xds,
		float *yds,
		float *zds,
		float imgXoffset,
		float imgYoffset,
		float imgZoffset,
		float *viewangles,
		float *zshifts,
		int nrviews,
		float *sinogram,
		int nrcols,         // image
		int nrrows,         //    does NOT
		int nrplanes,       //        contain a dummy 1 pixel frame
		void *pOrig,
		float vox_xy_size, //added fields
		float vox_z_size,//added fields 
		byte* xy_mask); //added fields for xy_mask

	void DD3Back_struct(float x0,
		float y0,
		float z0,
		int nrdetcols,
		int nrdetrows,
		float *xds,
		float *yds,
		float *zds,
		float imgXoffset,
		float imgYoffset,
		float imgZoffset,
		float *viewangles,
		float *zshifts,
		int nrviews,
		float *sinogram,
		int nrcols,         // image
		int nrrows,         //    does NOT
		int nrplanes,       //        contain a dummy 1 pixel frame
		void *pOrig,
		float vox_xy_size, //added fields
		float vox_z_size,//added fields 
		byte* xy_mask); //added fields for xy_mask
	//////////////////////////////////////////////////////////////////////////	
}

#endif

