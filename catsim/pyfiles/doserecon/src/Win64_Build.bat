ECHO OFF
REM ---------------------------------------------------------------------
REM Program Name: MakeWindows 
REM Purpose: Compliles C code into Windows run-time library (.dll)
REM Author:  Samit Basu Jed Pack Paul FitzGerald(GE Global Research)
REM Organization: GE Global Research and GE Healthcare, General Electric Company
REM Version:  6.0
REM Date:  Dec 30, 2014
REM Class: GE Confidential. General Electric Proprietary Data (c) 2012 General Electric Company 
REM ---------------------------------------------------------------------
ECHO ON

set PATH=C:\mingw64\bin;

cd C:\xxx\src\
C:\mingw64\bin\mingw32-make -f MakeWindows64 clean 
C:\mingw64\bin\mingw32-make -f MakeWindows64

@ECHO Windows build complete
@PAUSE
