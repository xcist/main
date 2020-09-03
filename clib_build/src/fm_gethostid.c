// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include <stdlib.h>
#include <unistd.h>

/*
 * gethostid() returns a long, but the id is a 32-bit int (the man page does not specify if it is
 * signed or unsigned).  I type convert the long to an unsigned int so that I can return it to
 * FreeMat (FreeMat does not support return types of int64 or uint64 as of this writing).  The
 * IQTB code apparently converts the result of gethostid() to an unsigned long.  If this turns
 * out to be a problem, we can change the call semantics so that FreeMat passes a uint64 by reference
 * to this function and this function populates it.
 */
unsigned int fm_gethostid(void)
{
    if (sizeof(long) == 8) {
      return (unsigned int) (0x00000000FFFFFFFF & gethostid());
    } else {
      return (unsigned int) gethostid();
    }
}


