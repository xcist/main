
import numpy as np
from stl import mesh

def Phantom_CAD_to_Polygonal(cadFilename, phantomFilename=None, materialName='water',
                             materialId=1, mesReductionRatio=-1, sizeScale=0,
                             offcenter=(0, 0, 0), drawCAD=True):

    if phantomFilename is None:
        dotId = cadFilename.rfind('.')
        phantomFilename = cadFilename[:dotId] + '.ppm'

    if mesReductionRatio > 0:
        # Load the STL CAD file
        meshStl = mesh.Mesh.from_file(cadFilename)
        
        # Reduce facet number
        numReducedFaces = int(meshStl.vectors.shape[0] * mesReductionRatio)
        meshStl.reduce(numReducedFaces)

        # Get the vertices and faces
        v = np.array(meshStl.vectors.reshape(-1, 3), dtype=np.float32)
        f = np.arange(v.shape[0]).reshape(-1, 3)

    else:
        # Load the STL CAD file without reduction
        v, f, _, _ = mesh.read_surface(cadFilename)

    # Off-center
    if any(offcenter):
        v += np.array(offcenter, dtype=np.float32)

    # Rescale
    if sizeScale != 0:
        v *= sizeScale

    # Show CAD
    if drawCAD:
        from mpl_toolkits import mplot3d
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = mplot3d.Axes3D(fig)
        ax.add_collection3d(mplot3d.art3d.Poly3DCollection(v[f], facecolors=(0.7, 0.7, 0.75), edgecolors='k', linewidths=0.2))
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.view_init(elev=20, azim=30)
        plt.show()

    # Save polygon phantom
    materialList = ['water', 'plexi', 'al', 'ti', 'fe', 'cu', 'graphite']
    materialList[materialId] = materialName

    with open(phantomFilename, 'w') as fid:
        fid.write('%% Mesh Phantom for polygonal projector\n')
        fid.write('%% Author: Mingye Wu, GE GRC\n\n')
        fid.write('obj=[];object=[];                                             %% DO NOT EDIT THIS LINE\n')
        fid.write('update = ''[obj,object]=AddObject(obj,object,materialList);'';  %% DO NOT EDIT THIS LINE\n\n')

        fid.write('materialList = { ')
        for mat in materialList:
            fid.write("'%s' " % mat)
        fid.write('};\n\n')

        fid.write('obj.density = 1;\n')
        fid.write('obj.material = %d;\n' % materialId)

        fid.write('obj.vertices = [')
        for vert in v:
            fid.write('%f %f %f;' % tuple(vert))
        fid.write('];\n')

        fid.write('obj.tri_inds = [')
        for face in f:
            fid.write('%d %d %d;' % tuple(face))
        fid.write('];\n')

        fid.write('eval(update);\n')

    return phantomFilename

