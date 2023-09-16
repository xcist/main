#ifndef P_NLOG_INLINE_H
#define P_NLOG_INLINE_H
// Note: You must compile this file with low optimization level (-O0 or -O1 on gcc 3.4.6)
// or the consts defined below (e.g., coeff1, etc) will be optimized out because they are
// only reference in the asm block.  There is probably a way around this (other than a
// gratuitous use of the variables), but I'm sure not how.
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#ifdef __APPLE__
const unsigned int coeff1[] __attribute__((aligned(16))) = { 0x3f7fffc4,0x3f7fffc4,0x3f7fffc4,0x3f7fffc4 };
const unsigned int coeff2[] __attribute__((aligned(16))) = { 0xbeffef80,0xbeffef80,0xbeffef80,0xbeffef80 };
const unsigned int coeff3[] __attribute__((aligned(16))) = { 0x3ea9e190,0x3ea9e190,0x3ea9e190,0x3ea9e190 };
const unsigned int coeff4[] __attribute__((aligned(16))) = { 0xbe7682ec,0xbe7682ec,0xbe7682ec,0xbe7682ec };
const unsigned int coeff5[] __attribute__((aligned(16))) = { 0x3e2bad82,0x3e2bad82,0x3e2bad82,0x3e2bad82 };
const unsigned int coeff6[] __attribute__((aligned(16))) = { 0xbdc33c0e,0xbdc33c0e,0xbdc33c0e,0xbdc33c0e };
const unsigned int coeff7[] __attribute__((aligned(16))) = { 0x3d13d187,0x3d13d187,0x3d13d187,0x3d13d187 };
const unsigned int coeff8[] __attribute__((aligned(16))) = { 0xbbd37841,0xbbd37841,0xbbd37841,0xbbd37841 };
const unsigned int or_mask[] __attribute__((aligned(16))) = { 0x3f800000,0x3f800000,0x3f800000,0x3f800000 };
const unsigned int and_mask[] __attribute__((aligned(16))) = { 0x807fffff,0x807fffff,0x807fffff,0x807fffff };
const unsigned int exp_sub[] __attribute__((aligned(16))) = { 127,127,127,127 };
const unsigned int exp_mult[] __attribute__((aligned(16))) = { 0x3f317218,0x3f317218,0x3f317218,0x3f317218 };
const unsigned int exp_extract[] __attribute__((aligned(16))) = { 0x7f800000,0x7f800000,0x7f800000,0x7f800000 };
#else
static const unsigned int coeff1[] __attribute__((aligned(16))) = { 0x3f7fffc4,0x3f7fffc4,0x3f7fffc4,0x3f7fffc4 };
static const unsigned int coeff2[] __attribute__((aligned(16))) = { 0xbeffef80,0xbeffef80,0xbeffef80,0xbeffef80 };
static const unsigned int coeff3[] __attribute__((aligned(16))) = { 0x3ea9e190,0x3ea9e190,0x3ea9e190,0x3ea9e190 };
static const unsigned int coeff4[] __attribute__((aligned(16))) = { 0xbe7682ec,0xbe7682ec,0xbe7682ec,0xbe7682ec };
static const unsigned int coeff5[] __attribute__((aligned(16))) = { 0x3e2bad82,0x3e2bad82,0x3e2bad82,0x3e2bad82 };
static const unsigned int coeff6[] __attribute__((aligned(16))) = { 0xbdc33c0e,0xbdc33c0e,0xbdc33c0e,0xbdc33c0e };
static const unsigned int coeff7[] __attribute__((aligned(16))) = { 0x3d13d187,0x3d13d187,0x3d13d187,0x3d13d187 };
static const unsigned int coeff8[] __attribute__((aligned(16))) = { 0xbbd37841,0xbbd37841,0xbbd37841,0xbbd37841 };
static const unsigned int or_mask[] __attribute__((aligned(16))) = { 0x3f800000,0x3f800000,0x3f800000,0x3f800000 };
static const unsigned int and_mask[] __attribute__((aligned(16))) = { 0x807fffff,0x807fffff,0x807fffff,0x807fffff };
static const unsigned int exp_sub[] __attribute__((aligned(16))) = { 127,127,127,127 };
static const unsigned int exp_mult[] __attribute__((aligned(16))) = { 0x3f317218,0x3f317218,0x3f317218,0x3f317218 };
static const unsigned int exp_extract[] __attribute__((aligned(16))) = { 0x7f800000,0x7f800000,0x7f800000,0x7f800000 };
#endif

void p_nlog_inline(float *dest, float *src, int len);
#endif // P_NLOG_INLINE_H
