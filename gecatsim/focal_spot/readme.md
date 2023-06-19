This folder contains focal spot images in `.npz` format for XCIST. Currently it only has VCT focal spot images.

The `.npz` format should have following keys:
* 'data': 2d focal spot images
* 'pixsize_x': pixel size in mm in the X direction
* 'pixsize_z': pixel size in mm in the Z direction
* [optional]


Example usages of the new focal spot model are:
1. Uniform:
`scanner.focalspotCallback = "SetFocalspot"
scanner.focalspotShape = "Uniform"
scanner.focalspotWidth = 1.0
scanner.focalspotLength = 1.0`
2. Gaussian
`scanner.focalspotCallback = "SetFocalspot"
scanner.focalspotShape = "Gaussian"
scanner.focalspotWidth = 1.0
scanner.focalspotLength = 1.0
scanner.focalspotWidthThreshold = 0.5
scanner.focalspotLengthThreshold = 0.5`
3. Customized focal spot image
`scanner.focalspotCallback = "SetFocalspot"
scanner.focalspotData = "xxx/xxx/xxx.npz"
scanner.focalspotWidth = 1.0
scanner.focalspotLength = 1.0
scanner.focalspotPixSizeX = 0.04
scanner.focalspotPixSizeZ = 0.04
#scanner.focalspotWidthThreshold = 0.5
scanner.focalspotLengthThreshold = 0.5`
