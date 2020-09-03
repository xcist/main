// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>

/*
 * This code is written to match the IQTB hiq_uid executable in terms of returning a PID.
 * hiq_uid returns the PID of the process, which is basically a random number.  We can't
 * just use getpid() in FreeMat because it will always return the PID of FreeMat.  Thus,
 * we fork off a process just to get its PID, which is a bit strange.
 */
pid_t fm_getrandpid(void)
{
    pid_t rc, child;

    if ( (child = fork()) != 0 ) {
        if ( (rc = wait(NULL)) != child) {
            fprintf(stderr, "Error: wait returned %d instead of %d in fm_getrandpid.  Was that a zombie process?\r\n", rc, child);
            exit(-1);
        }
        return child;
    } else {
        /* I am the child.  My only purpose was to generate a PID, which is a bit sad. */
        /* Do _not_ call just exit() because we don't want to call things registered for execution at atexit(). */
        _exit(0);
    }        
}


