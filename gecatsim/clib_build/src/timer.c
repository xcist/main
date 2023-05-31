// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#undef WIN32

#ifdef WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

#ifdef WIN32
#define ULONGLONG unsigned __int64
#define DWORD     unsigned
#else
#define ULONGLONG unsigned long long
#endif

#ifdef WIN32
ULONGLONG __inline cycle_counter(void)
{
  union
  {
    ULONGLONG		u64;
    struct
    {
      DWORD		low;
      DWORD		high;
    }				dw;
  }					cycles;
  
  __asm
    {
      _emit	0x0f;			
      _emit	0x31;
      mov	dword ptr cycles.dw.low, eax;
      mov	dword ptr cycles.dw.high,edx;
    }
  return cycles.u64;					
}

#else
ULONGLONG __inline cycle_counter(void)
{
  unsigned long low, high;
  __asm__ __volatile__("rdtsc;"
		       "movl %%eax, %0;"
		       "movl %%edx, %1"
		       : "=r" (low), "=r" (high)
		       : /* no input */
		       : "eax", "edx");
  return (((ULONGLONG) high << 32) | (ULONGLONG) low);
}
#endif

#ifdef WIN32
__declspec(dllexport)
#endif
double x86timer() {
  ULONGLONG count;
  count = cycle_counter();
  return (double) count;
}


