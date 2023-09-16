// Note: You must compile this file with low optimization level (-O0 or -O1 on gcc 3.4.6)
// or the consts defined below (e.g., coeff1, etc) will be optimized out because they are
// only reference in the asm block.  There is probably a way around this (other than a
// gratuitous use of the variables), but I'm sure not how.
#include "p_nlog_inline.h"
#ifdef WIN32
#define posix_memalign(p, a, s) (((*(p)) = _aligned_malloc((s), (a))), *(p) ?0 :errno)
#endif

void p_nlog_inline(float *dest, float *src, int len)
{
    // Make sure that len is divisible by four.  We could handle the non divisible by
    // four case by using the scalar SSE instruction equivalents, but for now we don't
    // worry about it.
    assert( (len & 0x3) == 0 );

    // We use unaligned loads and stores, so it is not necessary to enforce alignment
    // on the input and output.  It should be more efficient to use aligned loads and
    // stores, so if efficiency is required, then aligned can be required and the
    // two movups instructions can be changed to movaps.
    // We use the rip to generate position independent code on x86_64.
#ifndef __APPLE__
    __asm__( ".intel_syntax noprefix\n\t"
#elif defined __APPLE__
    __asm__( 
#endif
         "sar          %0,2\n\t"
#ifdef __i386__
         "movaps       xmm3,[or_mask]\n\t"
         "movaps       xmm4,[coeff6]\n\t"
         "movaps       xmm5,[coeff5]\n\t"
         "movaps       xmm6,[coeff4]\n\t"
         "movaps       xmm7,[coeff2]\n\t"
         "repeat:\n\t"
         "movups       xmm1,[%1]\n\t"
         "movaps       xmm0,xmm1\n\t"
         "andps        xmm1,[and_mask]\n\t"
         "orps         xmm1,xmm3\n\t"
         "movaps       xmm2,[coeff8]\n\t"
         "subps        xmm1,xmm3\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,[coeff7]\n\t"
         "psrld        xmm0,23\n\t"
         "mulps        xmm2,xmm1\n\t"
         "psubd        xmm0,[exp_sub]\n\t"
         "addps        xmm2,xmm4\n\t"
         "cvtdq2ps     xmm0,xmm0\n\t"
         "mulps        xmm2,xmm1\n\t"
         "mulps        xmm0,[exp_mult]\n\t"
         "addps        xmm2,xmm5\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,xmm6\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,[coeff3]\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,xmm7\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,[coeff1]\n\t"
#elif __x86_64__  &&  !defined __APPLE__
         "movaps       xmm3,[rip+or_mask]\n\t"
         "movaps       xmm4,[rip+coeff6]\n\t"
         "movaps       xmm5,[rip+coeff5]\n\t"
         "movaps       xmm6,[rip+coeff4]\n\t"
         "movaps       xmm7,[rip+coeff2]\n\t"
         "repeat:\n\t"
         "movups       xmm1,[%1]\n\t"
         "movaps       xmm0,xmm1\n\t"
         "andps        xmm1,[rip+and_mask]\n\t"
         "orps         xmm1,xmm3\n\t"
         "movaps       xmm2,[rip+coeff8]\n\t"
         "subps        xmm1,xmm3\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,[rip+coeff7]\n\t"
         "psrld        xmm0,23\n\t"
         "mulps        xmm2,xmm1\n\t"
         "psubd        xmm0,[rip+exp_sub]\n\t"
         "addps        xmm2,xmm4\n\t"
         "cvtdq2ps     xmm0,xmm0\n\t"
         "mulps        xmm2,xmm1\n\t"
         "mulps        xmm0,[rip+exp_mult]\n\t"
         "addps        xmm2,xmm5\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,xmm6\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,[rip+coeff3]\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,xmm7\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,[rip+coeff1]\n\t"
#elif __x86_64__ && defined __APPLE__
         "movaps       xmm3,[rip+_or_mask]\n\t"
         "movaps       xmm4,[rip+_coeff6]\n\t"
         "movaps       xmm5,[rip+_coeff5]\n\t"
         "movaps       xmm6,[rip+_coeff4]\n\t"
         "movaps       xmm7,[rip+_coeff2]\n\t"
         "repeat:\n\t"
         "movups       xmm1,[%1]\n\t"
         "movaps       xmm0,xmm1\n\t"
         "andps        xmm1,[rip+_and_mask]\n\t"
         "orps         xmm1,xmm3\n\t"
         "movaps       xmm2,[rip+_coeff8]\n\t"
         "subps        xmm1,xmm3\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,[rip+_coeff7]\n\t"
         "psrld        xmm0,23\n\t"
         "mulps        xmm2,xmm1\n\t"
         "psubd        xmm0,[rip+_exp_sub]\n\t"
         "addps        xmm2,xmm4\n\t"
         "cvtdq2ps     xmm0,xmm0\n\t"
         "mulps        xmm2,xmm1\n\t"
         "mulps        xmm0,[rip+_exp_mult]\n\t"
         "addps        xmm2,xmm5\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,xmm6\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,[rip+_coeff3]\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,xmm7\n\t"
         "mulps        xmm2,xmm1\n\t"
         "addps        xmm2,[rip+_coeff1]\n\t"
#endif
         "mulps        xmm2,xmm1\n\t"
         "add          %1,16\n\t"
         "addps        xmm2,xmm0\n\t"
         "subps        xmm0,xmm0\n\t"
         "subps        xmm0,xmm2\n\t"
         "movups       [%2],xmm0\n\t"
         "add          %2,16\n\t"
         "dec          %0\n\t"
         "jnz          repeat\n\t"
         ".att_syntax\n\t"
         : /* no output registers */
         : "r" (len), "r" (src), "r" (dest) : "memory" /* %0 == len, %1 == src, %2 == dest */  );
}

// The following is some code to test the p_nlog_inline function versus p_nlog.  It
// should be commented out unless testing is needed.  The cycle_counter functions
// depend on the code that reads the hardware cycle counter using rdtsc -- I haven't
// checked this code in yet because I haven't tested it on many architectures.  Anyway,
// it's not needed by the main code.

#if 0

#define NELEM (1024*1024*16)

int main(int argc, char **argv)
{
    float *f, *g1, *g2;
    unsigned long long t1, t2;
    double t;
    int i;

    if (posix_memalign( (void *) &f, 16, NELEM*sizeof(float)) != 0) {
        fprintf(stderr, "Error: posix_memalign failed on line %d in file %s\r\n", __LINE__, __FILE__);
        exit(-1);
    }
    if (posix_memalign( (void *) &g1, 16, NELEM*sizeof(float)) != 0) {
        fprintf(stderr, "Error: posix_memalign failed on line %d in file %s\r\n", __LINE__, __FILE__);
        exit(-1);
    }
    if (posix_memalign( (void *) &g2, 16, NELEM*sizeof(float)) != 0) {
        fprintf(stderr, "Error: posix_memalign failed on line %d in file %s\r\n", __LINE__, __FILE__);
        exit(-1);
    }

    for (i = 0; i < NELEM; i++) {
        f[i] = 1 + drand48() * 100;
    }

    for (i = 0; i < 10; i++) {
        memset(g1, 0, NELEM*sizeof(float));
        t1 = cycle_counter();
        p_nlog_inline(g1, f, NELEM);
        t2 = cycle_counter();
        printf("Neglog took %llu cycles.\r\n", t2-t1);
    }

    for (i = 0; i < 10; i++) {
        memset(g2, 0, NELEM*sizeof(float));
        t1 = cycle_counter();
        p_nlog(g2, f, NELEM);
        t2 = cycle_counter();
        printf("Orig neglog took %llu cycles.\r\n", t2-t1);
    }

    for (i = 0; i < NELEM; i++) {
        if (g1[i] != g2[i]) {
            printf("Error: g1[%d] = %f != g2[%d] = %f\r\n", i, g1[i], i, g2[i]);
           exit(-1);
        }
    }

    return 0;
}

#endif
