// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* global storage for large vectors */
#define MAX_FFT_SIZE 2*2048*16
#define POLY_ORDER 4

/*
**      TWIDDLE FACTOR TABLES FOR NEW FFT ROUTINES
*/
double twiddle[3*MAX_FFT_SIZE/4];
int bit_rev = 1;

double freqWindow[MAX_FFT_SIZE];
double starterKernel[MAX_FFT_SIZE];
double *starterKernelFft;
double reconFilterFft[MAX_FFT_SIZE];
double viewFft[MAX_FFT_SIZE];
double *tempArr;

void make_par_starter (int interp , int fftsiz , float detsp , double *startk) {
  double sampleSpacing;
  int hfftsz, iloop;
  
  /* Clear starter kernel */
  for ( iloop = 0; iloop < fftsiz; iloop++ )
    startk[iloop] = 0.0;
  
  hfftsz = fftsiz / 2;
  sampleSpacing = detsp / (double) interp;
  
  
  startk[0] = M_PI / (8.0 * sampleSpacing);
  for ( iloop = 1; iloop <= hfftsz; iloop += 2 ) {
    startk[iloop] = -1.0 / (2.0 * M_PI * sampleSpacing * iloop * iloop);
    startk[fftsiz - iloop] = startk[iloop];
  }
    
  for(iloop=0;iloop<fftsiz;iloop++)  {
    /* startk[iloop] = 2.0 * (double) interp * sampleSpacing * startk[iloop]; */
    startk[iloop] = 2.0 * (double) 2.0 * sampleSpacing * startk[iloop];
  }
  return;
}


#ifdef M_PI
#define TWO_PI 2.0*M_PI
#else
#define TWO_PI 6.2831853071795862
#endif

#ifdef M_SQRT1_2
#define RSQRT2 M_SQRT1_2
#else
#define RSQRT2 0.70710678118654752440
#endif


/*
**	rvfft.c -- Real valued, in-place, Cooley-Tukey 
**			radix-2 FFT routine
**
**	Real input and output in array x
**
**	Length is N = 2 ** m
**
**	Decimation in time, cos/sin in innermost loop
**
**	Output in order:
**	[Re(0),Re(1), ... ,Re(N/2),Im(N/2-1), ... , Im(1)]
**
**	The sine/cosine table for the twiddle factors is
**	expected to be supplied in the following format:
**
**	twiddle[0]       = sin(0*2*pi/n)
**	twiddle[1]       = sin(1*2*pi/n)
**	...
**	twiddle[n/2 - 1] = sin((n/2-1)*2*pi/n)
**
**	This corresponds to the first half period of a sine wave.
**
**	Based on Sorensen et al.: Real Valued FFT Algorithms
**	IEEE Transactions on Acoustics, Speech, and Signal Processing
**	Vol. ASSP-35, No. 6, June 1987
*/

void rvfft( double x[], int n, int m, int bitrev, int twiddleNSize, double twiddle[] ){

  int i,j,k;
  int n1,n2,n4;
  int i1, i2, i3, i4;
  int stage;
  int skip;
  int separation;
  int cosStart;
  int sinStart;
  int twiddleRatio;
  double xt;
  double cc,ss;
  double t1,t2;

  /*
  **	CALCULATE THE RATIO OF TWIDDLE FACTOR TABLE
  */
  twiddleRatio = twiddleNSize/n;

  /*
  **	separation between the sin and cos tables
  **	is the transform length divided by 4
  */
  separation = (n >> 2) * twiddleRatio;
  /*
  **	CHECK TO SEE IF WE DO THE BIT REVERSAL
  */
  n1 = n-1;
  n2 = (n >> 1);

  if (bitrev)
    {
      j = 0;
	
      for (i=1; i < n1; i++)
	{
	  k = n2;
	  while (k <=j )
	    {
	      j = j - k;
	      k = (k >> 1);
	    }

	  j += k;

	  if (i < j)
	    {
	      xt   = x[j];
	      x[j] = x[i];
	      x[i] = xt;
	    }
	}
    }


  /*
  **	PERFORM THE LENGTH 2 BUTTERFLIES
  **	STAGE 1
  */

  for (i=0; i<n; i+=2)
    {
      xt     = x[i];
      x[i]   = xt + x[i+1];
      x[i+1] = xt - x[i+1];
    }

  /*
  **	PERFORM THE REST OF THE STAGES
  **	STAGE 2 -> STAGE M
  */

  n2 = 1;

  for (stage = 2; stage <= m; stage++)
    {

      n4   = n2;
      n2   = (n4 << 1);
      n1   = (n2 << 1);
      skip = (n >> stage) * twiddleRatio;


      for (i=0; i<n; i+=n1)
	{
	  xt         = x[i];
	  x[i]       = xt + x[i+n2];
	  x[i+n2]    = xt - x[i+n2];
	  x[i+n4+n2] = -x[i+n4+n2];

	  sinStart = skip;
	  cosStart = skip + separation;
	  /*
	  **	NOTE: THE FIRST RUN THROUGH N4 = 1
	  */
	  for (j = 1; j < n4; j++)
	    {
	      i1 = i + j;
	      i2 = i - j + n2;
	      i3 = i + j + n2;
	      i4 = i - j + n1;
				
	      cc = twiddle[cosStart];
	      ss = twiddle[sinStart];

	      sinStart += skip;
	      cosStart += skip;

	      t1 = x[i3]*cc + x[i4]*ss;
	      t2 = x[i3]*ss - x[i4]*cc;

	      x[i4] =  x[i2] - t2;
	      x[i3] = -x[i2] - t2;
	      x[i2] =  x[i1] - t1;
	      x[i1] =  x[i1] + t1;
	    }
	}
    }
}

/*
**	irvfft.c -- Real valued, in-place, split-radix
**			IFFT routine
**
**	Hermitian symmetric imput and real output in array x
**
**	Length is n = 2 ** m
**
**	Decimation in frequency, cos/sin in second loop
**
**	Input order:
**	[Re(0),Re(1), ... ,Re(N/2),Im(N/2-1), ... , Im(1)]
**
**	The sine/cosine table for the twiddle factors are
**	expected to be supplied in the following order:
**
**      twiddle[0]       = sin(0*2*pi/n)
**      twiddle[1]       = sin(1*2*pi/n)
**      ...
**      twiddle[3*n/4 - 1] = sin((3*n/4-1)*2*pi/n)
**
**      This corresponds to the first 3/4 period of a sine wave.
**
**      Based on Sorensen et al.: Real Valued FFT Algorithms
**      IEEE Transactions on Acoustics, Speech, and Signal Processing
**      Vol. ASSP-35, No. 6, June 1987
*/

void irvfft(double x[], int n, int m, int bitrev, double twiddle[]) {

  int i,j,k;
  int n2,n4,n8;
  int i0,i1,i2,i3,i4,i5,i6,i7,i8;
  int is,id;
  int n1;

  int sinStart, sin3Start;
  int cosStart, cos3Start;
  int startAngle;
  int separation;
	
  double ss1, ss3;
  double cc1, cc3;
  double t1,t2,t3,t4,t5;
  double xt;
  double r1;
  double scale;

  n2 = n << 1;

  separation = n/4;


  for (k = 1; k < m; k++)
    {
      is = 0;
      id = n2;
      n2 = n2 >> 1;
      n4 = n2 >> 2;
      n8 = n4 >> 1;

      startAngle = n/n2;

      do
	{
	  for (i=is; i<n; i+=id)
	    {
	      i1 = i;
	      i2 = i1 + n4;
	      i3 = i2 + n4;
	      i4 = i3 + n4;

	      t1    = x[i1] - x[i3];
	      x[i1]  = x[i1] + x[i3];
	      x[i2] =         2.0 * x[i2];
	      x[i3] = t1    - 2.0 * x[i4];
	      x[i4] = t1    + 2.0 * x[i4];

	      if (n4 != 1)
		{
		  i1 += n8;
		  i2 += n8;
		  i3 += n8;
		  i4 += n8;

		  t1 = (x[i2] - x[i1]) * RSQRT2;
		  t2 = (x[i4] + x[i3]) * RSQRT2;

		  x[i1] = x[i1] + x[i2];
		  x[i2] = x[i4] - x[i3];
		  x[i3] = 2.0 * (-t2-t1);
		  x[i4] = 2.0 * (-t2+t1);
		}
	    }

	  is = (id << 1) - n2;
	  id = (id << 2);

	} while ( is < (n-1) );

      sinStart = startAngle;
      cosStart = sinStart + separation;



      for (j=2; j <= n8; j++)
	{

	  sin3Start = 3*sinStart;
	  cos3Start = sin3Start + separation;

	  cc1 = twiddle[cosStart];
	  ss1 = twiddle[sinStart];
	  cc3 = twiddle[cos3Start];
	  ss3 = twiddle[sin3Start];

	  sinStart = j * startAngle;
	  cosStart = sinStart + separation;

	  is = 0;

	  id = 2 * n2;
	  do
	    {
	      for (i=is; i<n; i+=id)
		{
		  i1 = i + j - 1;
		  i2 = i1 + n4;
		  i3 = i2 + n4;
		  i4 = i3 + n4;
		  i5 = i + n4 - j + 1;
		  i6 = i5 + n4;
		  i7 = i6 + n4;
		  i8 = i7 + n4;

		  t1    = x[i1] - x[i6];
		  x[i1] = x[i1] + x[i6];

		  t2    = x[i5] - x[i2];
		  x[i5] = x[i2] + x[i5];

		  t3    = x[i8] + x[i3];
		  x[i6] = x[i8] - x[i3];

		  t4    = x[i4] + x[i7];
		  x[i2] = x[i4] - x[i7];

		  t5    = t1 - t4;
		  t1    = t1 + t4;
		  t4    = t2 - t3;
		  t2    = t2 + t3;

		  x[i3] =  t5 * cc1 + t4 * ss1;
		  x[i7] = -t4 * cc1 + t5 * ss1;
		  x[i4] =  t1 * cc3 - t2 * ss3;
		  x[i8] =  t2 * cc3 + t1 * ss3;
		}

	      is = 2 * id - n2;
	      id = 4 * id;

	    } while ( is < (n-1) );

	}
    }

  /*
  **	LENGTH TWO BUTTERFLIES
  */

  is = 0;
  id = 4;
	
  do
    {

      for ( i0 = is; i0 < n; i0+= id)
	{
	  i1 = i0 + 1;
	  r1 = x[i0];

	  x[i0] = r1 + x[i1];
	  x[i1] = r1 - x[i1];
	}

      is = (id << 1) - 2;
      id = (id << 2);

    } while ( is < n );

  /*
  **      CHECK TO SEE IF WE DO THE BIT REVERSAL
  */
  n1 = n-1;
  n2 = (n >> 1); 

  if (bitrev)
    {
      j = 0;

      for (i=1; i < n1; i++)
	{
	  k = n2;
	  while (k <=j )
	    {
	      j = j - k;
	      k = (k >> 1);
	    }
 
	  j += k;
 
	  if (i < j)
	    {
	      xt   = x[j];
	      x[j] = x[i];
	      x[i] = xt;
	    }
	}
    }

  /*
  **	SCALE BY 1/N
  */
  scale = 1.0/n;

  for (i=0; i<n; i++)
    x[i] *= scale;

}
 
 
/**********************************
 **	Calculate twiddle factors
 **********************************/

void initrealroots(double *roots, unsigned int size) {
  int k;

  *roots++ = 1.0;

  for (k=1; k< 3*size/4; k++)
    {
      *roots++ = sin(TWO_PI/(double)size * (double)k);
    }
}

void make_starter(int interp, int fftsiz, double detsp, double stoi, double startk[])
{
  double alpha, numerator, denominator;
  int hfftsz, iloop;

  /* Clear starter kernel */
  for ( iloop = 0; iloop < fftsiz; iloop++ )
    startk[iloop] = 0.0;

  hfftsz = fftsiz / 2;
  alpha  = detsp / (double) interp / stoi;
  numerator = -2.0 * alpha * alpha / M_PI;
  /*
    printf ("----- half fft size = %d\n", hfftsz);
    printf ("----- alpha = %e  and numerator = %e\n", alpha, numerator);
  */
  startk[0] = M_PI / 2.0;
  for ( iloop = 1; iloop <= hfftsz; iloop += 2 )
    {
      denominator = sin( (double) iloop * alpha );
      denominator *= denominator;
      startk[iloop] = numerator / denominator;
      startk[fftsiz - iloop] = startk[iloop];
      /*	printf ("----- iloop = %d   startk = %e\n", iloop, startk[iloop]); */
    }
}


void make_window (float fcutof[], int porder, float mpolyc[], int wintyp, int interp, int fftsiz, double mwindw[]) {
  double rbase, rtemp;
  int halfft, nramp, nstart, nend, iloop, polyloop;

  /* Determine start and end array indices for basic window array */
  halfft = fftsiz/2;
  nstart = (int) ( fcutof[0] * (float) (halfft + 1) + 0.9999999);
  if (nstart < 0) nstart = 0;
  nend = (int) ( fcutof[1] * (float) (halfft + 1) + 0.9999999);
  if (nend > halfft) nend = halfft;
  nramp = nend - nstart;

  /* Clear basic window array */
  for ( iloop=0; iloop <fftsiz; iloop++)
    mwindw[iloop] = 0.0;

  /* Fill basic window array with ones up to and including the
     element following the start index */
  for ( iloop = 0; iloop <= nstart; iloop++)
    mwindw[iloop] = 1.0;

  rbase = M_PI / (double) nramp;
  /*
    printf ("** Value of window step increment = %f\n", rbase);
  */

  /* If window type == SINC, calculate a sync basic window */
  if (wintyp == 1)
    for ( iloop = 1; iloop < nramp; iloop++ )
      {
	rtemp = rbase * (float) iloop;
	mwindw[nstart + iloop] = sin(rtemp) / rtemp;
      }

  /* Else assume window type == HANNING */
  else
    for ( iloop = 0; iloop < nramp; iloop++ )
      mwindw[nstart + iloop] = 0.5*(1.0 + cos(rbase * (float) iloop));

  rbase = 2.0 / (double) (fftsiz / interp);

  /* Calculate the polynomial window */
  for ( iloop = 0; iloop < halfft; iloop++ )
    {
      rtemp = mpolyc[0];
      for ( polyloop = 1; polyloop <= porder; polyloop++ )
	rtemp = rtemp*rbase*iloop +  mpolyc[polyloop];
      mwindw[iloop] *= (double)(interp * interp) * rtemp;
    }

  /* copy left half of window to righ half */
  for ( iloop = 1; iloop < halfft; iloop++ )
    mwindw[halfft+iloop] = mwindw[halfft-iloop];

}

#ifdef WIN32
__declspec(dllexport)
#endif
void gemsrampcurved(float*sino_in, float *sino_out, int windowType, 
		    int freqInterp, float dfov, float detSpacing, 
		    float fpixf, float fmax, float flow, float polyCoefs[5],
		    int numRows, int numCols, int numViews, float SID, int isParallel) {
  	int numPtsToConvolve = numCols;
  	int numDets = numPtsToConvolve;
  	int ifftSize;
  	/* find smallest power of 2 >= to numDets */
  	for (ifftSize = 1; ifftSize<numDets; ifftSize = ifftSize * 2);
  	ifftSize = ifftSize * 2;	/* float fft size to avoid aliasing
				   due to circular convolution */
  	ifftSize = ifftSize * freqInterp;
  	/* account for zero padding */
  	/* cutoff frequency calculations */
  	float sampleSizeAtIso = detSpacing;
  	float pixelSize = dfov / 512; /* ASSUME 512 image size */
  	float pixelCutoff = sampleSizeAtIso / pixelSize;
  	float cutoffs[2];
  	/* high cutoff */
  	if (fmax < fpixf*pixelCutoff)
    	cutoffs[1] = fmax;
  	else
    	cutoffs[1] = fpixf*pixelCutoff;
  	/* low cutoff */
  	cutoffs[0] = flow * cutoffs[1];
  	/* divide cutoff freq's by freqInterp factor to avoid aliasing
     due to the high freq we may have injected by adding
     zeros between points */
  	cutoffs[1] = cutoffs[1] / freqInterp;
  	cutoffs[0] = cutoffs[0] / freqInterp;
  	int fftSize = ifftSize / freqInterp;
  	int numConvolvedPts = numPtsToConvolve * freqInterp;

	if (isParallel == 0) {
		make_starter(freqInterp, ifftSize, detSpacing, SID, starterKernel);
	} else {
		make_par_starter(freqInterp, ifftSize, sampleSizeAtIso, starterKernel);
	}
  	/*
  	**      NEW FFT
  	**      INIT REAL ROOTS FOR FFT ROUTINES
  	*/
  	initrealroots(twiddle, ifftSize);
  
  	int ifftSizeD2 = ifftSize/2;
  	int fftSizeD2 = fftSize/2;

  	/* generate frequency domain "window" function */
  	make_window(cutoffs, POLY_ORDER, polyCoefs, windowType, freqInterp, ifftSize, freqWindow);
  
  	int log2ifftSize = (int)( log((double)ifftSize) / log(2.0) + 0.0001 );
  	int log2fftSize = log2ifftSize - freqInterp/2;

  	starterKernelFft = &starterKernel[0];
  
  	rvfft(starterKernelFft, ifftSize, log2ifftSize,
		bit_rev, ifftSize, twiddle);

  	/* multiply window and fft of starter kernel and scale
     (scale factor is non-unity only for parallel beam case)
     note that imaginary portion of fft should be 0
     since starterKernel is symmetric
  	*/
  	int i, j;
  	int viewNum, rows;
  	for (i=0; i<= ifftSize/2; i++)  {
    	reconFilterFft[i] = freqWindow[i] * starterKernelFft[i];
  	}

  	/*  what is ifft?   */

  	/***    'FILTERing VIEW LOOP'    ***/
  
  	/*  loop through projections.  viewNum is 0-based and relative
      to first view to process.
  	*/
  
  	for (viewNum=0; viewNum<numViews; viewNum++)
    {
      for(rows = 0; rows < numRows; rows++)
		{
	  	memset(viewFft, '\0', ifftSize*sizeof(double));
	  	for (i=0; i<numPtsToConvolve; i++)
	    	viewFft[i] = (double) sino_in[viewNum*numRows*numCols + rows*numCols + i];
	  	/*
	  	**      TAKE FORWARD FFT
	  	*/
	  	rvfft(viewFft, fftSize, log2fftSize, bit_rev, ifftSize, twiddle);
	  	/*      Since we really wanted to inject zeros between every point
		  in the projection, use the DSP property that when zeros
		  are injected, the effect is to replicate the data in freq.
		  space.  The fft would be periodic w/ freqInterp periods.
	  	*/
	  	/*      freq. pt. @ Nyquist is a special case when we do interpolation
		  since it was origninally scaled by fftSize (not by fftSize/2
		  like the other non-zero frequencies)
	  	*/
	  	if (freqInterp > 1)
	    	viewFft[fftSizeD2] *= 0.5;
	  	/*
	  	**      NOW REPLICATE THE SPECTRUM freqInterp-1 TIMES
	  	*/
	  	tempArr = &freqWindow[0];
	  	/*
	  	**      REPLICATE THE SPECTRUM
	  	*/
	  	if (freqInterp > 1) {
	    	for (j=0; j<freqInterp/2; j++) {                
	      	int tempindex,tempindex1, tempindex2, tempindex3;
	      	int index1, index2;
	      	tempindex = j*fftSize;
	      	tempindex1 = (j+1)*fftSize;
	      	tempindex2 = freqInterp/4+2;
	      	tempindex3 = freqInterp/4+1;
	      	for (i=1; i<fftSizeD2; i++) {
			/*
			**                      THIS IS THE REAL PART
			*/
			index1 = i + tempindex;
			index2 = i;
			tempArr[index1] = viewFft[index2];
			index1 = tempindex1 - i;
			tempArr[index1] = viewFft[index2];
			/*
			**                      THIS IS THE IMAGINARY PART
			*/
			index1 = tempindex+(tempindex2*fftSize)-i;
			index2 = fftSize - i;
			tempArr[index1] = viewFft[index2];
			index1 = tempindex+(tempindex3*fftSize)+i;
			tempArr[index1] = -viewFft[index2];
	      	}
	      	/*
	      	**              THESE ARE THE SPECIAL CASES FOR THE REAL PART
	      	*/
	      	tempArr[j*fftSize] = viewFft[0];
	      	tempArr[(j+1)*fftSize] = viewFft[0];
	      	tempArr[(j*fftSize)+fftSizeD2] = 2.0*viewFft[fftSizeD2];
	      	/*
	      	**              THESE ARE THE SPECIAL CASES FOR THE IMAGINARY PART
	      	*/
	      	tempArr[3*ifftSize/4] = 0.0;
	      	index1 = (unsigned int)((0.5* (freqInterp + 1) + j) *fftSize);
	      	tempArr[index1] = 0.0;
	      
	    	}
	  	}
	  	else    /* no interpolation */
	    {
	    	memcpy(tempArr, viewFft, fftSize*sizeof(double));
	    }
	  	/*
	  	** MULTIPLY BY THE FILTER AND TAKE THE INVERSE TRNASFORM
	  	** WE TAKE INTO ACCOUNT THAT THE IMAGINARY PART OF THE
	  	** RECON FILTER IS ZERO.
	  	**      (a+bi)(c+di) = ac + bci  when d = 0
	  	*/
	  	viewFft[0] = reconFilterFft[0] * tempArr[0];
	  	for (i=1; i < ifftSizeD2; i++)
	    {
	      	viewFft[i] = reconFilterFft[i] * tempArr[i];
	    }
	  	viewFft[ifftSizeD2] = reconFilterFft[ifftSizeD2] * tempArr[ifftSizeD2];
	  	for (i=1; i < ifftSizeD2; i++)
	    {
	      	viewFft[ifftSizeD2+i] = reconFilterFft[ifftSizeD2-i] * tempArr[ifftSizeD2+i];
	    } 
	  	/*
	  	**  TAKE INVERSE FFT
	  	*/
	  	irvfft(viewFft, ifftSize, log2ifftSize, bit_rev, twiddle);
	  	/*
	  	**  CONVERT TO FLOAT FROM DOUBLE
	  	*/
	  	for (i=0; i<numConvolvedPts; i++)
	    	sino_out[viewNum*numRows*freqInterp*numCols + rows*numCols*freqInterp + i] = (float)(viewFft[i]);
		}
    }  /* end FILTERing VIEW LOOP */
}


#ifdef WIN32
__declspec(dllexport)
#endif
void gemsramp(float*sino_in, float *sino_out, int windowType, 
	      int freqInterp, float dfov, float detSpacing, 
	      float fpixf, float fmax, float flow, float polyCoefs[5],
	      int numRows, int numCols, int numViews) {
  	int numPtsToConvolve = numCols;
  	int numDets = numPtsToConvolve;
  	int ifftSize;
  	/* find smallest power of 2 >= to numDets */
  	for (ifftSize = 1; ifftSize<numDets; ifftSize = ifftSize * 2);
  	ifftSize = ifftSize * 2;	/* float fft size to avoid aliasing
				   due to circular convolution */
  	ifftSize = ifftSize * freqInterp;
  	/* account for zero padding */
  	/* cutoff frequency calculations */
  	float sampleSizeAtIso = detSpacing;
  	float pixelSize = dfov / 512; /* ASSUME 512 image size */
  	float pixelCutoff = sampleSizeAtIso / pixelSize;
  	float cutoffs[2];
  	/* high cutoff */
  	if (fmax < fpixf*pixelCutoff)
    	cutoffs[1] = fmax;
  	else
    	cutoffs[1] = fpixf*pixelCutoff;
  	/* low cutoff */
  	cutoffs[0] = flow * cutoffs[1];
  	/* divide cutoff freq's by freqInterp factor to avoid aliasing
     due to the high freq we may have injected by adding
     zeros between points */
  	cutoffs[1] = cutoffs[1] / freqInterp;
  	cutoffs[0] = cutoffs[0] / freqInterp;
  	int fftSize = ifftSize / freqInterp;
  	int numConvolvedPts = numPtsToConvolve * freqInterp;

	make_par_starter(freqInterp, ifftSize, detSpacing, starterKernel);
  	/*
  	**      NEW FFT
  	**      INIT REAL ROOTS FOR FFT ROUTINES
  	*/
  	initrealroots(twiddle, ifftSize);
  
  	int ifftSizeD2 = ifftSize/2;
  	int fftSizeD2 = fftSize/2;

  	/* generate frequency domain "window" function */
  	make_window(cutoffs, POLY_ORDER, polyCoefs, windowType, freqInterp, ifftSize, freqWindow);

  	int log2ifftSize = (int)( log((double)ifftSize) / log(2.0) + 0.0001 );
  	int log2fftSize = log2ifftSize - freqInterp/2;

  	starterKernelFft = &starterKernel[0];
  
  	rvfft(starterKernelFft, ifftSize, log2ifftSize,
		bit_rev, ifftSize, twiddle);

  	/* multiply window and fft of starter kernel and scale
     (scale factor is non-unity only for parallel beam case)
     note that imaginary portion of fft should be 0
     since starterKernel is symmetric
  	*/
  	int i, j;
  	int viewNum, rows;
  	for (i=0; i<= ifftSize/2; i++)  {
    	reconFilterFft[i] = freqWindow[i] * starterKernelFft[i];
 	}

  	/*  what is ifft?   */

  	/***    'FILTERing VIEW LOOP'    ***/
  
  	/*  loop through projections.  viewNum is 0-based and relative
      to first view to process.
  	*/
  
  	for (viewNum=0; viewNum<numViews; viewNum++)
    {
      for(rows = 0; rows < numRows; rows++)
		{
	  	memset(viewFft, '\0', ifftSize*sizeof(double));
	  	for (i=0; i<numPtsToConvolve; i++)
	    	viewFft[i] = (double) sino_in[viewNum*numRows*numCols + rows*numCols + i];
	  	/*
	  	**      TAKE FORWARD FFT
	  	*/
	  	rvfft(viewFft, fftSize, log2fftSize, bit_rev, ifftSize, twiddle);
	  	/*      Since we really wanted to inject zeros between every point
		  in the projection, use the DSP property that when zeros
		  are injected, the effect is to replicate the data in freq.
		  space.  The fft would be periodic w/ freqInterp periods.
	  	*/
	  	/*      freq. pt. @ Nyquist is a special case when we do interpolation
		  since it was origninally scaled by fftSize (not by fftSize/2
		  like the other non-zero frequencies)
	  	*/
	  	viewFft[fftSizeD2] *= 0.5;
	  	/*
	  	**      NOW REPLICATE THE SPECTRUM freqInterp-1 TIMES
	  	*/
	  	tempArr = &freqWindow[0];
	  	/*
	  	**      REPLICATE THE SPECTRUM
	  	*/
	  	if (freqInterp > 1) {
	    	for (j=0; j<freqInterp/2; j++) {                
	      	int tempindex,tempindex1, tempindex2, tempindex3;
	      	int index1, index2;
	      	tempindex = j*fftSize;
	      	tempindex1 = (j+1)*fftSize;
	      	tempindex2 = freqInterp/4+2;
	      	tempindex3 = freqInterp/4+1;
	      	for (i=1; i<fftSizeD2; i++) {
			/*
			**                      THIS IS THE REAL PART
			*/
			index1 = i + tempindex;
			index2 = i;
			tempArr[index1] = viewFft[index2];
			index1 = tempindex1 - i;
			tempArr[index1] = viewFft[index2];
			/*
			**                      THIS IS THE IMAGINARY PART
			*/
			index1 = tempindex+(tempindex2*fftSize)-i;
			index2 = fftSize - i;
			tempArr[index1] = viewFft[index2];
			index1 = tempindex+(tempindex3*fftSize)+i;
			tempArr[index1] = -viewFft[index2];
	      	}
	      	/*
	      	**              THESE ARE THE SPECIAL CASES FOR THE REAL PART
	      	*/
	      	tempArr[j*fftSize] = viewFft[0];
	      	tempArr[(j+1)*fftSize] = viewFft[0];
	      	tempArr[(j*fftSize)+fftSizeD2] = 2.0*viewFft[fftSizeD2];
	      	/*
	      	**              THESE ARE THE SPECIAL CASES FOR THE IMAGINARY PART
	      	*/
	      	tempArr[3*ifftSize/4] = 0.0;
	      	index1 = (unsigned int)((0.5* (freqInterp + 1) + j) *fftSize);
	      	tempArr[index1] = 0.0;
	      
	    	}
	  	}
	  	else    /* no interpolation */
	    {
	      	memcpy(tempArr, viewFft, fftSize*sizeof(double));
	    }
	  	/*
	  	** MULTIPLY BY THE FILTER AND TAKE THE INVERSE TRNASFORM
	  	** WE TAKE INTO ACCOUNT THAT THE IMAGINARY PART OF THE
	  	** RECON FILTER IS ZERO.
	  	**      (a+bi)(c+di) = ac + bci  when d = 0
	  	*/
	  	viewFft[0] = reconFilterFft[0] * tempArr[0];
	  	for (i=1; i < ifftSizeD2; i++)
	    {
	      	viewFft[i] = reconFilterFft[i] * tempArr[i];
	    }
	  	viewFft[ifftSizeD2] = reconFilterFft[ifftSizeD2] * tempArr[ifftSizeD2];
	  	for (i=1; i < ifftSizeD2; i++)
	    {
	      	viewFft[ifftSizeD2+i] = reconFilterFft[ifftSizeD2-i] * tempArr[ifftSizeD2+i];
	    } 
	  	/*
	  	**  TAKE INVERSE FFT
	  	*/
	  	irvfft(viewFft, ifftSize, log2ifftSize, bit_rev, twiddle);
	  	/*
	  	**  CONVERT TO FLOAT FROM DOUBLE
	  	*/
	  	for (i=0; i<numConvolvedPts; i++)
	    	sino_out[viewNum*numRows*freqInterp*numCols + rows*numCols*freqInterp + i] = (float)(viewFft[i]);
		}
    }  /* end FILTERing VIEW LOOP */
}

/*
 * Includes the erroneous PI frequency
 */
void gemsramp_old(float*sino_in, float *sino_out, int windowType, 
		  int freqInterp, float dfov, float detSpacing, 
		  float fpixf, float fmax, float flow, float polyCoefs[5],
		  int numRows, int numCols, int numViews) {

  	int numPtsToConvolve = numCols;
  	int numDets = numPtsToConvolve;
  	int ifftSize;
  	/* find smallest power of 2 >= to numDets */
  	for (ifftSize = 1; ifftSize<numDets; ifftSize = ifftSize * 2);
  	ifftSize = ifftSize * 2;	/* float fft size to avoid aliasing
				   due to circular convolution */
  	ifftSize = ifftSize * freqInterp;
  	/* account for zero padding */
  	/* cutoff frequency calculations */
  	float sampleSizeAtIso = detSpacing;
  	float pixelSize = dfov / 512; /* ASSUME 512 image size */
  	float pixelCutoff = sampleSizeAtIso / pixelSize;
  	float cutoffs[2];
  	/* high cutoff */
  	if (fmax < fpixf*pixelCutoff)
    	cutoffs[1] = fmax;
  	else
    	cutoffs[1] = fpixf*pixelCutoff;
  	/* low cutoff */
  	cutoffs[0] = flow * cutoffs[1];
  	/* divide cutoff freq's by freqInterp factor to avoid aliasing
     due to the high freq we may have injected by adding
     zeros between points */
  	cutoffs[1] = cutoffs[1] / freqInterp;
  	cutoffs[0] = cutoffs[0] / freqInterp;
  	int fftSize = ifftSize / freqInterp;
  	int numConvolvedPts = numPtsToConvolve * freqInterp;

  	make_par_starter(freqInterp, ifftSize, detSpacing, starterKernel);
  	/*
  	**      NEW FFT
  	**      INIT REAL ROOTS FOR FFT ROUTINES
  	*/
  	initrealroots(twiddle, ifftSize);
  
  	int ifftSizeD2 = ifftSize/2;
  	int fftSizeD2 = fftSize/2;

  	/* generate frequency domain "window" function */
  	make_window(cutoffs, POLY_ORDER, polyCoefs, windowType, freqInterp, ifftSize, freqWindow);

  	int log2ifftSize = (int)( log((double)ifftSize) / log(2.0) + 0.0001 );
  	int log2fftSize = log2ifftSize - freqInterp/2;

  	starterKernelFft = &starterKernel[0];
  
  	rvfft(starterKernelFft, ifftSize, log2ifftSize,
		bit_rev, ifftSize, twiddle);

  	/* multiply window and fft of starter kernel and scale
     (scale factor is non-unity only for parallel beam case)
     note that imaginary portion of fft should be 0
     since starterKernel is symmetric
  	*/
  	int i, j;
  	int viewNum, rows;
  	for (i=0; i<= ifftSize/2; i++)  {
    	reconFilterFft[i] = freqWindow[i] * starterKernelFft[i];
  	}

  	/*  what is ifft?   */

  	/***    'FILTERing VIEW LOOP'    ***/
  
  	/*  loop through projections.  viewNum is 0-based and relative
      to first view to process.
  	*/
  
  	for (viewNum=0; viewNum<numViews; viewNum++)
    {
      for(rows = 0; rows < numRows; rows++)
		{
	  	memset(viewFft, '\0', ifftSize*sizeof(double));
	  	for (i=0; i<numPtsToConvolve; i++)
	    	viewFft[i] = (double) sino_in[viewNum*numRows*numCols + rows*numCols + i];
	  	/*
	  	**      TAKE FORWARD FFT
	  	*/
	  	rvfft(viewFft, fftSize, log2fftSize, bit_rev, ifftSize, twiddle);
	  	/*      Since we really wanted to inject zeros between every point
		  in the projection, use the DSP property that when zeros
		  are injected, the effect is to replicate the data in freq.
		  space.  The fft would be periodic w/ freqInterp periods.
	  	*/
	  	/*      freq. pt. @ Nyquist is a special case when we do interpolation
		  since it was origninally scaled by fftSize (not by fftSize/2
		  like the other non-zero frequencies)
	  	*/
	  	viewFft[fftSizeD2] *= 0.5;
	  	/*
	  	**      NOW REPLICATE THE SPECTRUM freqInterp-1 TIMES
	  	*/
	  	tempArr = &freqWindow[0];
	  	/*
	  	**      REPLICATE THE SPECTRUM
	  	*/
	  	if (freqInterp > 1) {
	    	for (j=0; j<freqInterp/2; j++) {                
	      	int tempindex,tempindex1, tempindex2, tempindex3;
	      	int index1, index2;
	      	tempindex = j*fftSize;
	      	tempindex1 = (j+1)*fftSize;
	      	tempindex2 = freqInterp/4+2;
	      	tempindex3 = freqInterp/4+1;
	      	for (i=1; i<fftSizeD2; i++) {
			/*
			**                      THIS IS THE REAL PART
			*/
			index1 = i + tempindex;
			index2 = i;
			tempArr[index1] = viewFft[index2];
			index1 = tempindex1 - i;
			tempArr[index1] = viewFft[index2];
			/*
			**                      THIS IS THE IMAGINARY PART
			*/
			index1 = tempindex+(tempindex2*fftSize)-i;
			index2 = fftSize - i;
			tempArr[index1] = viewFft[index2];
			index1 = tempindex+(tempindex3*fftSize)+i;
			tempArr[index1] = -viewFft[index2];
	      	}
	      	/*
	      	**              THESE ARE THE SPECIAL CASES FOR THE REAL PART
	      	*/
	      	tempArr[j*fftSize] = viewFft[0];
	      	tempArr[(j+1)*fftSize] = viewFft[0];
	      	tempArr[(j*fftSize)+fftSizeD2] = viewFft[fftSizeD2];
	      	/*
	      	**              THESE ARE THE SPECIAL CASES FOR THE IMAGINARY PART
	      	*/
	      	tempArr[3*ifftSize/4] = 0.0;
	      	index1 = (unsigned int)((0.5* (freqInterp + 1) + j) *fftSize);
	      	tempArr[index1] = 0.0;
	      
	    	}
	  	}
	  	else    /* no interpolation */
	    {
	      	memcpy(tempArr, viewFft, fftSize*sizeof(double));
	    }
	  	/*
	  	** MULTIPLY BY THE FILTER AND TAKE THE INVERSE TRNASFORM
	  	** WE TAKE INTO ACCOUNT THAT THE IMAGINARY PART OF THE
	  	** RECON FILTER IS ZERO.
	  	**      (a+bi)(c+di) = ac + bci  when d = 0
	  	*/
	  	viewFft[0] = reconFilterFft[0] * tempArr[0];
	  	for (i=1; i < ifftSizeD2; i++)
	    {
	      	viewFft[i] = reconFilterFft[i] * tempArr[i];
	    }
	  	viewFft[ifftSizeD2] = reconFilterFft[ifftSizeD2] * tempArr[ifftSizeD2];
	  	for (i=1; i < ifftSizeD2; i++)
	    {
	      	viewFft[ifftSizeD2+i] = reconFilterFft[ifftSizeD2-i] * tempArr[ifftSizeD2+i];
	    } 
	  	/*
	  	**  TAKE INVERSE FFT
	  	*/
	  	irvfft(viewFft, ifftSize, log2ifftSize, bit_rev, twiddle);
	  	/*
	  	**  CONVERT TO FLOAT FROM DOUBLE
	  	*/
	  	for (i=0; i<numConvolvedPts; i++)
	    	sino_out[viewNum*numRows*freqInterp*numCols + rows*numCols*freqInterp + i] = (float)(viewFft[i]);
		}
    }  /* end FILTERing VIEW LOOP */
}


