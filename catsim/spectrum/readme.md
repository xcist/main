# How to generate a customied spectrum
1. use `XCISTspectrum.m` to generate spectrum at a user specified kVp/angle/energybin size. Note we will have a `.mat` file that is not readable to XCIST
2. modify the spectrum filename in `convert.py` and run it to obtain a spectrum file ending in `.dat` that can be read by XCIST
