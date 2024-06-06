
# How to compile recon lib

# For Linux
Go to the `reconstruction\src` folder and run following command:

`gcc -fPIC -fopenmp -shared -o fdk_equiAngle.so interface_fdk_angle.c`

This will generate a `fdk_equiAngle.so` file and copy it to your `reconstruction\lib` folder.

# For Mac OS
Same as for Linux but make sure to use the "gcc" compiler and not the "clang" compiler.

You may need to install gcc manually (e.g., brew install gcc).
If needed, specify the full path to the gcc compiler when compiling, e.g.:

`/usr/local/bin/gcc-14 -fPIC -fopenmp -shared -o fdk_equiAngle.so interface_fdk_angle.c`

# For Windows
in cmd, run `xxx\mingw64\bin\gcc.exe -Wall -O -g -fPIC -shared -o fdk_equiAngle.dll .\interface_fdk_angle.c `
