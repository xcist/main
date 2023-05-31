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

"C:\Program Files\CodeBlocks\MinGW\bin"\mingw32-make -f MakeWindows32 clean
"C:\Program Files\CodeBlocks\MinGW\bin"\mingw32-make -f MakeWindows32

@ECHO Windows build complete
@PAUSE
