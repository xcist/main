// Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

//#include "negative_log.h"
#include <cmath>
#include <string.h>
#include <stdlib.h>
#ifdef WIN32
#define DLLEXPORT __declspec(dllexport)
//https://stackoverflow.com/questions/33696092/whats-the-correct-replacement-for-posix-memalign-in-window
#define posix_memalign(p, a, s) (((*(p)) = _aligned_malloc((s), (a))), *(p) ?0 :errno)
#else
#define DLLEXPORT
#endif
//#include "p_nlog_inline.h"
#define USE_C_LOG

#define NLOG_ZERO_REPLACEMENT     1.175494351E-38F
static void negative_log_with_das_underrange_corr(int row_count, int col_count, float *iview, float *oview);

static void negative_log_without_das_underrange_corr(int row_count, int col_count, float *input, float *output);

static void log_replace_lte_zero(float *input, float zero_replacement, int length);

extern "C" void p_nlog_inline(float*, float*, int);
void min_vector(float *input, float *min, int length) {
  float local_min = input[0];
  for (int i = 1; i < length; i++)
    if (input[i] < local_min)
      local_min = input[i];
  *min = local_min;
}

void nlog(float *dest, float *src, int numpoints) {
  #ifdef USE_C_LOG
    for (int i = 0; i < numpoints; i++)
      dest[i] = -logf(src[i]);
  #else
    static void *source_aligned, *dest_aligned;
    static int aligned_size = 0;
    // it would be better to cache these arrays and only allocate if/when the size changes
    if (aligned_size != numpoints) {
      if (aligned_size != 0) {
        free(source_aligned);
        free(dest_aligned);
      }
      posix_memalign(&source_aligned, 32, numpoints*sizeof(float));
      posix_memalign(&dest_aligned, 32, numpoints*sizeof(float));
      aligned_size = numpoints;
    }
    memcpy(source_aligned, src, numpoints*sizeof(float));
    p_nlog_inline((float *)dest_aligned, (float *)source_aligned, numpoints);
    memcpy(dest, dest_aligned, numpoints*sizeof(float));
  #endif
}

// ***********************************************************************
// Function applyAAR_AARi (5-point kernel)
// ***********************************************************************
int applyAAR_AARi(float *input, float *scratch, float lowClip, int length) {
  float alpha;
  float alphaInv;
  float fiveAlpha;
  float noiseChan;
  float smoothingGain;
  float x;
  float viewError;
  int count = 0;
  int lenm2;
  int i;

  // INITIALIZE
  alpha = 0.001f;
  alphaInv = 1.0f/alpha;
  fiveAlpha = 5.0f * alpha;
  lenm2 = length - 2;

  // COPY INPUT INTO SCATCH BUFFER
  for (i = 0; i < length; i++)
    scratch[i] = input[i];

  // SMOOTH OUT LOW SIGNAL AREAS
  for (i = 2; i < lenm2; i++) {
    if (input[i] <= fiveAlpha) {
      // HIGH PASS FILTER USING 5 POINT KERNEL
      noiseChan = 0.625f*scratch[i] - 0.25f*(scratch[i+1] + scratch[i-1]) - 0.0625f*(scratch[i+2] + scratch[i-2]);

      // SMOOTHING GAIN IS FIFTH ORDER POLYNOMIAL ESTIMATE OF exp()
      x = input[i] * alphaInv;
      smoothingGain = ((((-0.000973364f*x + 0.01664038f)*x - 0.1181265f)*x + 0.4528541f)*x - 0.9823636f)*x + 0.9990775f;
      if (smoothingGain > 1.0f)
        smoothingGain = 1.0f;

      // COMPUTE VIEW ERROR
      viewError = noiseChan * smoothingGain;

      // SUBTRACT ERROR FROM SIGNAL
      input[i] -= viewError;
      count++;
    }
  }
  return(count);
}


// ***********************************************************************
// Function applyAAR_AARi2 (7-point kernel)
// This is a newer version of improved AAR for some improvement in the IQ
// ***********************************************************************
int applyAAR_AARi2(float *input, float *scratch, float lowClip, int length) {
  float alpha;
  float alphaInv;
  float fiveAlpha;
  float noiseChan;
  float smoothingGain;
  float x;
  float viewError;
  int count = 0;
  int lenm3;
  int i;

  // INITIALIZE
  alpha = 0.001f;
  alphaInv = 1.0f/alpha;
  fiveAlpha = 5.0f * alpha;
  lenm3 = length - 3;

  // COPY INPUT INTO SCATCH BUFFER
  for (i = 0; i < length; i++)
    scratch[i] = input[i];

  // SMOOTH OUT LOW SIGNAL AREAS
  for (i = 3; i < lenm3; i++) {
    if (input[i] <= fiveAlpha) {
      // HIGH PASS FILTER USING 7 POINT KERNEL
      noiseChan = scratch[i] - ((scratch[i] + scratch[i+1] + scratch[i-1] +
                  scratch[i+2] + scratch[i-2] + scratch[i+3] + scratch[i-3]) / 7.0f);

      // SMOOTHING GAIN IS FIFTH ORDER POLYNOMIAL ESTIMATE OF exp()
      x = input[i] * alphaInv;
      smoothingGain = ((((-0.000973364f*x + 0.01664038f)*x - 0.1181265f)*x + 0.4528541f)*x - 0.9823636f)*x + 0.9990775f;
      if (smoothingGain > 1.0f)
        smoothingGain = 1.0f;

      // COMPUTE VIEW ERROR
      viewError = noiseChan * smoothingGain;

      // SUBTRACT ERROR FROM SIGNAL
      input[i] -= viewError;
      count++;
    }
  }
  return(count);
}

// If das_underrange_corr is non-zero, then we do das underrange correction (AARi2).  Otherwise, we
// only replace pre-log values of <= 0 with NLOG_ZERO_REPLACEMENT
extern "C" {
DLLEXPORT
void negative_log(int row_count, int col_count, float *iview, float *oview, int das_underrange_corr) {
  if (das_underrange_corr) {
    negative_log_with_das_underrange_corr(row_count, col_count, iview, oview);
  } else {
    negative_log_without_das_underrange_corr(row_count, col_count, iview, oview);
  }
}
}


static void negative_log_with_das_underrange_corr(int row_count, int col_count, float *iview, float *oview) {
  float lowClip = expf(-13.0f);
  float lowThr = expf(-10.3f);
  float pThr = expf(-8.0f);
  int length = col_count;

  for (int row = 0; row < row_count; row++) {
    float *input = iview + row * col_count;
    float *output = oview + row * col_count;
    float min;
    int i;

    min_vector(input, &min, length);

    // Clip 3 edge channels, to be consistent with the product
    for (i = 0; i < 3; i++) {
      if (input[i] < lowClip)
        input[i] = lowClip;
      if (input[length-i-1] < lowClip)
        input[length-i-1] = lowClip;
    }

    if (min < lowThr) {
      for (i = 0; i < length; i++)
        if (input[i] <= lowClip)
          input[i] = lowClip * (1.0 + input[i] * 50.0f);
      applyAAR_AARi2(input, output, lowClip, length);
      min_vector(input, &min, length);
    }

    if (min < pThr)
      applyAAR_AARi(input, output, lowClip, length);

    for (i = 0; i < length; i++)
      if (input[i] < lowClip)
        input[i] = lowClip;

    nlog(output, input, length);
  }
}


static void negative_log_without_das_underrange_corr(int row_count, int col_count, float *input, float *output) {
  int length = row_count * col_count;
  float min;

  min_vector(input, &min, length);
  if (min <= 0)
    log_replace_lte_zero(input, NLOG_ZERO_REPLACEMENT, length);

  nlog(output, input, length);
}


// lte == less than or equal
static void log_replace_lte_zero(float *input, float zero_replacement, int length) {
  int i;

  for (i = 0; i < length; i++) {
    if (input[i] == 0)
      input[i] = zero_replacement;
    else if (input[i] < 0)
      input[i] = -input[i];
  }
}
