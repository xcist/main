
# How to compile in linux

Go to the `reconstruction\src` folder and run following command:

`gcc -fPIC -fopenmp -shared -o fdk_equiAngle.so interface_fdk_angle.c`

This will generate a `fdk_equiAngle.so` file and copy it to your 
