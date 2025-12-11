# How to compile on Windows

1. unzip `pthreads-2-9-1-release.zip` under `catsim\lib\`
2. modify lines 10 and 11 in `MakeWindoes64`:
```
LFLAGS = -LC:\mingw64\lib -LC:\Users\xxx\software\main-master\catsim\lib -lpthreadGC2_x64
IFLAGS = -IC:\mingw64\include -IC:\Users\xxx\software\main-master\catsim\lib\Pre-built.2\include
```
3. compile using `BuildWin64.bat`

# How to compile on Mac

Due to an issue when running multi-threaded version of ncat projector compiled with `gcc` on mac, we need to first install gcc, then replace `gcc` in `MakeVariable_3` ***and*** `/usr/local/opt/gcc/bin/gcc-14` in MakeMacOS with the path  to real gcc, such as `/usr/local/opt/gcc/bin/gcc-14` and then recompile. Note that the default `gcc` on mac may simply be an alias to `clang`. To verify this, run 'gcc --version' from the command line.

To compile, cd into the src directory and then run:

`make -f ../MakeMacOS`

