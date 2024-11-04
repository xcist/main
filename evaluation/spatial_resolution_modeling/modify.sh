#! /bin/bash

# gain factor
sed -i "s/scanner.detectionGain.*/scanner.detectionGain = 10/g" */*/*/*cfg
sed -i "s/scanner.detectionGain.*/scanner.detectionGain = 10/g" */*/*/*/*cfg
# enoise
sed -i "s/scanner.eNoise.*/scanner.eNoise = 1000/g" */*/*/*cfg
sed -i "s/scanner.eNoise.*/scanner.eNoise = 1000/g" */*/*/*/*cfg
