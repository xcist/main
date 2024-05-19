# How to compile on Windows

1. unzip `pthreads-2-9-1-release.zip` under `catsim\lib\`
2. modify lines 10 and 11 in `MakeWindoes64`:
```
LFLAGS = -LC:\mingw64\lib -LC:\Users\xxx\software\main-master\catsim\lib\Pre-built.2\dll\x64 -lpthreadGC2
IFLAGS = -IC:\mingw64\include -IC:\Users\xxx\software\main-master\catsim\lib\Pre-built.2\include
```
3. compile using `BuildWin64.bat`

# How to compile on Mac

Due to an issue when running multi-threaded version of ncat projector compiled with `gcc` on mac, we need to first install gcc, then replace `gcc` in `MakeVariable_3` with the path to path to real gcc, such as `/usr/local/opt/gcc/bin/gcc-14` and then recompile.b Note that the default `gcc` on mac is simply an alias to `clang`.
