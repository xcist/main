# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

def ReadMaterialFile(mtFile):
    # read material file
    d0 = []
    for line in open(mtFile, 'r', encoding='UTF-8'):
        line = line.strip() # remove the end '\n'
        if line and line[0].isdigit():
            d0.append(line)

    numberOfElements = int(d0[0])
    density = float(d0[1])
    atomicNumbers = []
    massFractions = []
    for ii in range(2, 2+numberOfElements):
        d1 = d0[ii].split()
        atomicNumbers.append(int(d1[0]))
        massFractions.append(float(d1[1]))
        
    # normalize mass fraction
    massSum = sum(massFractions)
    for ii in range(0, len(massFractions)):
        massFractions[ii] /= massSum

    return numberOfElements, density, atomicNumbers, massFractions

# if __name__ == '__main__':
#     (numberOfElements, density, atomicNumbers, massFractions) = ReadMaterialFile('water')
#     print(numberOfElements, density, atomicNumbers, massFractions)
