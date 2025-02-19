# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

"""
Aim
    This function converts a .pp phantom file to a .ppm file.

Inputs:
    PhantomFileBasename  The name of the phantom file, with no extension.
                         A file named [PhantomFileBasename '.pp'] must exist.
    Scale                (optional) A scale factor, to resize phantom, also used to convert from cm to mm (scale=10).
                         Default: 1

Outputs:
    A new file is written to [PhantomFileBasename '.ppm']
"""

import os
from math import *
from gecatsim.pyfiles.C_Phantom_Analytic_FORBILD_to_tmp import C_Phantom_Analytic_FORBILD_to_tmp

def phantom_analytic_pp_to_ppm(cfg, phantom_file_basename, scale=1):
    if not phantom_file_basename:
        raise ValueError('PhantomFileBasename must be specified')

    print(f'Converting phantom file {phantom_file_basename}.pp to .ppm file...')

    pp_phantom_filename = f'{phantom_file_basename}.pp'
    tmp_phantom_filename = f'{phantom_file_basename}.tmp'

    C_Phantom_Analytic_FORBILD_to_tmp(cfg, scale, pp_phantom_filename, tmp_phantom_filename)

    print(f'Reading and deleting {tmp_phantom_filename}.')
    lns = read_text_lines2(tmp_phantom_filename)
    os.remove(tmp_phantom_filename)
    ppm_phantom_filename = f'{phantom_file_basename}.ppm'
    print(f'Writing {ppm_phantom_filename}.')

    with open(ppm_phantom_filename, 'w') as fid:
        ind = [i for i, line in enumerate(lns) if line.startswith('#')]

        if len(ind) < 2:
            print('Warning: Material table not found. Using default material list.')
            fid.write('materialList = {\'DefaultMaterial\'};\n\n')
            offset = 0
        else:
            fid.write('materialList = {')

            for i in range(ind[1] - 1):
                mn = lns[i + 1].strip()
                fid.write(f"'{mn}' ")

            fid.write('};\n\n')
            offset = ind[1]
        i = 0

        while offset < len(lns) - 15:
            i += 1
            tmpo, tmpc, lns, offset = read_object(lns, offset)
            fid.write(f'object.center({i},:) = [{tmpo["cent"][0]} {tmpo["cent"][1]} {tmpo["cent"][2]}];\n')
            fid.write(f'object.half_axes({i},:) = [{tmpo["hax"][0]} {tmpo["hax"][1]} {tmpo["hax"][2]}];\n')
            fid.write(f'object.euler_angs({i},:) = [{tmpo["ea"][0]} {tmpo["ea"][1]} {tmpo["ea"][2]}];\n')
            fid.write(f'object.density({i}) = {tmpo["density"]};\n')
            fid.write(f'object.type({i}) = {tmpo["type"]};\n')
            fid.write(f'object.material({i}) = {tmpo["material"]};\n')
            fid.write(f'object.axial_lims({i},:) = [0 0];\n')
            fid.write(f'object.shape({i}) = 0;\n')
            fid.write(f'object.clip{{{i}}} = [')
            for j in range(len(tmpc)):
                fid.write(f'{tmpc[j][0]} {tmpc[j][1]} {tmpc[j][2]} {tmpc[j][3]};')
            fid.write('];\n\n')

    print(f'... done writing {ppm_phantom_filename}')

def read_text_lines2(file_name):
    try:
        with open(file_name, 'rb') as file:
            file_contents = file.read()
    except FileNotFoundError:
        raise FileNotFoundError(f"Error! Cannot find file {file_name}")

    line_indices = [0] + [i for i, byte in enumerate(file_contents) if byte == 10] + [len(file_contents)]
    max_line_width = max(line_indices[i+1] - line_indices[i] for i in range(len(line_indices) - 1))
    number_of_lines = len(line_indices) - 1
    lines = []

    for i in range(number_of_lines):
        line_length = line_indices[i+1] - line_indices[i] - 1
        if line_length > 0:
            line = file_contents[line_indices[i]:line_indices[i+1]].decode('utf-8').rstrip('\r')
            lines.append(line)

    return lines

def read_object(lines, offset):
    obj = {}
    clip = []

    obj['type'] = lines[0 + offset].split()[1]

    A = []
    for i in range(4):
        A.append([float(x) for x in lines[i + 1 + offset].split()])

    if sum([(A[3][i] - [0, 0, 0, 1][i]) ** 2 for i in range(4)]) > 1e-5:
        raise ValueError('Unexpected transform')

    obj['cent'] = [A[0][3], A[1][3], A[2][3]]
    B = [row[0:3] for row in A[0:3]]
    obj['hax'] = [sum([B[i][j] ** 2 for j in range(3)]) ** 0.5 for i in range(3)]
    B = [[B[i][j] / obj['hax'][i] for j in range(3)] for i in range(3)]
    B = list(map(list, zip(*B)))
    obj['ea'] = [angle * 180 / 3.141592653589793 for angle in euler_angs(B)]

    obj['transform'] = A
    tmp = lines[5 + offset].split()[0]

    if tmp == 'C':
        i = 1
        while tmp == 'C':
            v = [float(x) for x in lines[5 + offset][21:].split(',')]
            clip.append(v + [float(lines[6 + offset][26:])])
            clip[-1] = [-x for x in clip[-1]]
            offset += 3
            tmp = lines[5 + offset].split()[0]
            i += 1
    else:
        clip = []

    if obj['type'].startswith('Cyl'):
        obj['type'] = 2
        obj['hax'][2] /= 2
    elif obj['type'].startswith('Cub'):
        obj['type'] = 2
        Bt = list(map(list, zip(*B)))
        clip += [
            [Bt[j][i] * obj['hax'][i] / 2 + sum([obj['cent'][k] * Bt[j][k] for k in range(3)]) for i in range(3)] + [
                obj['hax'][i] / 2 - sum([obj['cent'][k] * Bt[j][k] for k in range(3)])] for j in range(2)]
        obj['hax'][:2] = [x * (2 ** 0.5) / 2 for x in obj['hax'][:2]]
    elif obj['type'].startswith('Sph'):
        obj['type'] = 1
    elif obj['type'].startswith('Ell'):
        obj_type = 1
    else:
        print('Unrecognized obj type')

    obj['material'] = int(lines[7 + offset].split()[1])
    obj['density'] = float(lines[8 + offset].split()[1])

    offset += 9

    return obj, clip, lines, offset

def euler_angs(R):
    p2 = acos(R[2][2])

    if abs(sin(p2)) > 1e-6:
        p1 = atan2(R[0][2] / sin(p2), R[1][2] / sin(p2))
        p3 = atan2(R[2][0] / sin(p2), -R[2][1] / sin(p2))
    elif cos(p2) > 0:
        p2, p3, p1 = 0, 0, atan2(R[0][1], R[0][0])
    else:
        p2, p3, p1 = pi, 0, atan2(-R[0][1], R[0][0])

    Rnew = [
        [cos(p1) * cos(p3) - sin(p1) * cos(p2) * sin(p3), cos(p1) * sin(p3) + sin(p1) * cos(p2) * cos(p3),
         sin(p1) * sin(p2)],
        [-sin(p1) * cos(p3) - cos(p1) * cos(p2) * sin(p3), -sin(p1) * sin(p3) + cos(p1) * cos(p2) * cos(p3),
         cos(p1) * sin(p2)],
        [sin(p2) * sin(p3), -sin(p2) * cos(p3), cos(p2)]
    ]

    err = sum([(Rnew[i][j] - R[i][j]) ** 2 for i in range(3) for j in range(3)]) ** 0.5

    if err > 1e-5:
        raise ValueError('Error computing Euler angles')

    return [p1, p2, p3]