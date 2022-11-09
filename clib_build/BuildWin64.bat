ECHO OFF
REM Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE
ECHO ON

REM show environment variables
set
CD clib_build\src

mingw32-make -f ..\MakeWindows64
if %ERRORLEVEL% NEQ 0 (
    exit /B %ERRORLEVEL%
)
move /Y libcatsim64.dll ..\..\catsim\lib
if %ERRORLEVEL% NEQ 0 (
    exit /B %ERRORLEVEL%
)

mingw32-make -f ..\MakeWindows64 clean
if %ERRORLEVEL% NEQ 0 (
    exit /B %ERRORLEVEL%
)

@ECHO Windows build complete
