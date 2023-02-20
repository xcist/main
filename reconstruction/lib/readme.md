
# How to compile recon lib

# for Linux
Go to the `reconstruction\src` folder and run following command:

`gcc -fPIC -fopenmp -shared -o fdk_equiAngle.so interface_fdk_angle.c`

This will generate a `fdk_equiAngle.so` file and copy it to your `reconstruction\lib` folder.

# For Windows
in cmd, run `xxx\mingw64\bin\gcc.exe -Wall -O -g -fPIC -shared -o fdk_equiAngle.dll .\interface_fdk_angle.c `
