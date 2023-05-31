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
