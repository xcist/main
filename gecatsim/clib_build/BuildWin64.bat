ECHO OFF
REM Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE
ECHO ON

set PATH=C:\mingw64\bin;%PATH%
CD src

C:\mingw64\bin\mingw32-make -f ..\MakeWindows64
move /Y libcatsim64.dll ..\..\lib
@PAUSE

C:\mingw64\bin\mingw32-make -f ..\MakeWindows64 clean

CD ..

@ECHO Windows build complete
@PAUSE
