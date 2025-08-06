// Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

#include <cmath>
#include <string.h>
#include <stdlib.h>
//#include "p_nlog_inline.h"

#define NLOG_ZERO_REPLACEMENT     1.175494351E-38F
void min_vector(float *input, float *min, int length);

void nlog(float *dest, float *src, int numpoints);

int applyAAR_AARi(float *input, float *scratch, float lowClip, int length);

int applyAAR_AARi2(float *input, float *scratch, float lowClip, int length);

void negative_log(int row_count, int col_count, float *iview, float *oview, int das_underrange_corr);

static void negative_log_with_das_underrange_corr(int row_count, int col_count, float *iview, float *oview);

static void negative_log_without_das_underrange_corr(int row_count, int col_count, float *input, float *output);

static void log_replace_lte_zero(float *input, float zero_replacement, int length);
