#ifndef __DD3__
#define __DD3__

// GE Proprietary

//extern "C"{

void DD3DoseRow(float imgX,
		float imgXstep,
		int nrcols,
		float imgZstart,
		float imgZstep,
		int nrplanes,
		float *pImg,
		float *pEnergymap, //NEW
		float *detX,
		int increment,
		float *detZ,
		float *scales,
		float z0,
		float *viewEnergyIn,
		float *viewEnergyDep,
		int nrdetrows,
		int nrdetcols,
		float *sils,
		float *pImgFirst,
		float *pImgLast,
		float *energyInFirst,
		float *energyInLast,
		float *energyDepFirst,
		float *energyDepLast,
		float *pEnergyMapFirst,
		float *pEnergyMapLast,
		float *detXFirst,
		float *detXLast,
		float *detZFirst,
		float *detZLast); //NEW
void DD3DoseView(float x0,
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
		    float *viewEnergyIn, // incoming energy : view with 1 pixel margin all around
		    float *viewEnergyDep, // depositied energy : empty view with 1 pixel margin all around
		    float *sils, // NEW slab-intersection-lengths - empty view with 1 pixel margin all around
		    int nrcols,
		    int nrrows,       // images
		    int nrplanes,     //     do NOT
		    float *pOrig,     //        contain
		    float *pTrans,//)    //             a dummy 1 pixel frame
		    float *pDose,//NEW
		    float *pDoseTrans,//NEW
		    float vox_xy_size,   // voxel size
		 float vox_z_size);    //field to represent vox sizes
void DD3Dose(float x0,
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
			float *pOrig,
		     float *pDose,// NEW
			float vox_xy_size, //added fields
		     float vox_z_size); //added fields 


void DD3FlatWBackRow(float imgX,
		float imgXstep,
		int nrcols,
		float imgZstart,
		float imgZstep,
		int nrplanes,
		float *pImg,
		float *detX,
		int increment,
		float *detZ,
		float *scales,
		float z0,
		float *view,
		 int nrdetrows);
void DD3FlatWBackView(float x0,
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
		 float dzdx,
		 float *sinogram,
		 float *view, // empty view with 1 pixel margin all around
		 int nrcols,
		 int nrrows,       // images
		 int nrplanes,     //     do NOT
		 float *pOrig,     //        contain
		 float *pTrans);    //             a dummy 1 pixel frame
void DD3FlatWBack(float x0,
	     float y0,
	     float z0,
	     int nrdetcols,
	     int nrdetrows,
	     float *xds,
	     float *yds,
	     float *zds,
	     float dzdx,
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
	      float *pOrig);
void DD3CurvedWBackRow(float imgX,
		float imgXstep,
		int nrcols,
		float imgZstart,
		float imgZstep,
		int nrplanes,
		float *pImg,
		float *detX,
		int increment,
		float *detZ,
		float *scales,
		float z0,
		float *view,
		 int nrdetrows);
void DD3CurvedWBackView(float x0,
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
		 float dzdx,
		 float *sinogram,
		 float *view, // empty view with 1 pixel margin all around
		 int nrcols,
		 int nrrows,       // images
		 int nrplanes,     //     do NOT
		 float *pOrig,     //        contain
		 float *pTrans);    //             a dummy 1 pixel frame
void DD3CurvedWBack(float x0,
	     float y0,
	     float z0,
	     int nrdetcols,
	     int nrdetrows,
	     float *xds,
	     float *yds,
	     float *zds,
	     float dzdx,
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
	      float *pOrig);


/*
 * DD3 transpose
 */
void DD3Transpose(int nrcols,
		  int nrrows,
		  int nrplanes,
		  float *pOrig,
		  float *pTrans);

/*
 * DD3 boundaries
 */
void DD3Boundaries(int nrBoundaries,
		   float *pCenters,
		   float *pBoundaries);

/*
 * DD3 add transpose
 */
void DD3AddTranspose(int nrcols,
		     int nrrows,
		     int nrplanes,
		     float *pOrig,
		     float *pTrans);

/*
 * DD3 row projector
 */
void DD3ProjRow(float imgX,
		float imgXstep,
		int nrcols,
		float imgZstart,
		float imgZstep,
		int nrplanes,
		float *pImg,
		float *detX,
		int increment,
		float *detZ,
		float *scales,
		float z0,
		float *view,
		int nrdetrows);

/*
 * DD3 row backprojector
 */
void DD3BackRow(float imgX,
		float imgXstep,
		int nrcols,
		float imgZstart,
		float imgZstep,
		int nrplanes,
		float *pImg,
		float *detX,
		int increment,
		float *detZ,
		float *scales,
		float z0,
		float *view,
		int nrdetrows);

/*
 * DD3 view projector
 */
void DD3ProjView(float x0,
		 float y0,
		 float z0,
		 int nrdetcols,
		 int nrdetrows,
		 int vertical,
		 float *xdi,
		 float *ydi,
		 float *detX,
		 float *detZ,
		 float *scales,
		 float dzdx,
		 float *view,
		 float *newproj,
		 int nrcols,
		 int nrrows,       // images
		 int nrplanes,     //     do NOT
		 float *pOrig,     //        contain
		 float *pTrans);   //             a dummy 1 pixel frame

/*
 * DD3 view backprojector
 */
void DD3BackView(float x0,
		 float y0,
		 float z0,
		 int nrdetcols,
		 int nrdetrows,
		 int vertical,
		 float *xdi,
		 float *ydi,
		 float *detX, // empty distance array (nrxdist + 2)
		 float *detZ,
		 float *scales,		
		 float dzdx,
		 float *sinogram,
		 float *view, // empty view with 1 pixel margin all around
		 int nrcols,
		 int nrrows,       // images
		 int nrplanes,     //     do NOT
		 float *pOrig,     //        contain
		 float *pTrans);   //             a dummy 1 pixel frame

/*
 * DD3 projector
 */

void DD3Proj(float x0,
	     float y0,
	     float z0,
	     int nrdetcols,
	     int nrdetrows,
	     float *xds,
	     float *yds,
	     float *zds,
	     float dzdx,
	     float imgXoffset,
	     float imgYoffset,
	     float imgZoffset,
	     float *viewangles,
	     float *zshifts,
	     int nrviews,
	     float *sinogram,
	     int nrcols,           // image
	     int nrrows,           //    does NOT
	     int nrplanes,         //        contain a dummy 1 pixel frame
	     float *pOrig);
void DD3Proj_mm(float x0,
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
	     int nrcols,           // image
	     int nrrows,           //    does NOT
	     int nrplanes,         //        contain a dummy 1 pixel frame
	     float *pOrig,
	     float voxel_xy_size,
	     float voxel_z_size);

/*
 * DD3 backprojector
 */
void DD3Back(float x0,
	     float y0,
	     float z0,
	     int nrdetcols,
	     int nrdetrows,
	     float *xds,
	     float *yds,
	     float *zds,
	     float dzdx,
	     float imgXoffset,
	     float imgYoffset,
	     float imgZoffset,
	     float *viewangles,
	     float *zshifts,
	     int nrviews,
	     float *sinogram,
	     int nrcols,           // image
	     int nrrows,           //    does NOT
	     int nrplanes,         //        contain a dummy 1 pixel frame
	     float *pOrig);
//}

#endif
