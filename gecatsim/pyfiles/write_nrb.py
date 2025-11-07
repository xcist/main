# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

def nrbwrite(filename, xp, yp, zp, uvecs, vvecs, mat_ind, obj_names):
    """
    Writes NURBS surface data to a file, including control points and knot vectors.

    Parameters:
    - filename (str): Output file name.
    - xp, yp, zp (list of np.ndarray): Control point coordinates.
    - uvecs, vvecs (list of np.ndarray): Knot vectors.
    - mat_ind (list of int): Material indices.
    - obj_names (list of str): Object names.
    """
    with open(filename, 'w') as fid:
        num_objs = len(xp)
        for i in range(num_objs):
            fid.write('\n')
            my_xp = xp[i].T
            my_yp = yp[i].T
            my_zp = zp[i].T
            uvec = uvecs[i]
            vvec = vvecs[i]
            Np, Nt = my_xp.shape
            fid.write(f'{obj_names[i]}\n')
            fid.write(f'{mat_ind[i]}\n')
            fid.write(f'{Np} :M\n')
            fid.write(f'{Nt} :N\n')
            fid.write('U Knot Vector\n')
            for val in uvec:
                fid.write(f'{val:1.6f}\n')
            fid.write('V Knot Vector\n')
            for val in vvec:
                fid.write(f'{val:1.6f}\n')
            fid.write('Control Points\n')
            for j in range(Np):
                for t in range(Nt):
                    fid.write(f'{my_xp[j, t]:1.6f} {my_yp[j, t]:1.6f} {my_zp[j, t]:1.6f}\n')
