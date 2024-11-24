# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
import numpy as np
import os
from gecatsim.pyfiles.C_Phantom_Analytic_FORBILD_to_tmp import C_Phantom_Analytic_FORBILD_to_tmp


def phantom_analytic_pp_to_ppm(phantom_file_basename, scale=1):

    if not phantom_file_basename:
        raise ValueError('PhantomFileBasename must be specified')

    pp_phantom_filename = f'{phantom_file_basename}.pp'
    tmp_phantom_filename = f'{phantom_file_basename}.tmp'

    C_Phantom_Analytic_FORBILD_to_tmp(scale, pp_phantom_filename, tmp_phantom_filename)

    lns = read_text_lines2(tmp_phantom_filename)
    os.remove(tmp_phantom_filename)

    ppm_phantom_filename = f'{phantom_file_basename}.ppm'

    with open(ppm_phantom_filename, 'w') as fid:
        ind = [i for i, line in enumerate(lns) if line.startswith('#')]
        if len(ind) < 2:
            raise ValueError('Material table not found')
        else:
            fid.write('materialList = {')
            for i in range(1, ind[1] - 1):
                mn = lns[i + 1].strip()
                fid.write(f"'{mn}' ")
            fid.write('};\n\n')
            offset = ind[1]

        i = 0
        while offset < len(lns) and lns[offset].strip():
            obj = {}
            obj['type'] = lns[10 + offset].strip()
            obj['name'] = lns[11 + offset].strip()
            obj['cent'] = np.array([float(x) for x in lns[12 + offset].strip().split()])
            obj['hax'] = np.array([float(x) for x in lns[13 + offset].strip().split()])
            A = np.array([float(x) for x in lns[14 + offset].strip().split()]).reshape(4, 4)

            if np.linalg.norm(A[3, :3]) > 1e-5:
                raise ValueError('Unexpected transform')

            obj['cent'] = A[:3, 3]
            B = A[:3, :3]
            obj['hax'] = np.sqrt(np.sum(B * B, axis=0))
            B = B / obj['hax']
            B = B.T
            obj['ea'] = euler_angs(B) * 180 / np.pi

            obj['transform'] = A
            tmp = lns[15 + offset].strip().split()[0]
            if tmp == 'C':  # clipping planes
                clip = []
                while tmp == 'C':
                    v = [float(x) for x in lns[15 + offset][21:].strip().split(',')]
                    clip.append(v + [float(lns[16 + offset][26:].strip())])
                    offset += 3
                    tmp = lns[15 + offset].strip().split()[0]
                clip = np.array(clip)
            else:
                clip = []

            if obj['type'].startswith('Cyl'):
                obj['type'] = 2
                obj['hax'][2] /= 2  # length needs to be divided by two for half axis
            elif obj['type'].startswith('Cub'):
                obj['type'] = 2  # change to cylinder
                Bt = B.T
                clip = np.vstack([clip, [Bt[:, 0], obj['hax'][0] / 2 + np.dot(obj['cent'], Bt[:, 0])]])
                clip = np.vstack([clip, [-Bt[:, 0], obj['hax'][0] / 2 - np.dot(obj['cent'], Bt[:, 0])]])
                clip = np.vstack([clip, [Bt[:, 1], obj['hax'][1] / 2 + np.dot(obj['cent'], Bt[:, 1])]])
                clip = np.vstack([clip, [-Bt[:, 1], obj['hax'][1] / 2 - np.dot(obj['cent'], Bt[:, 1])]])
                obj['hax'][:2] *= np.sqrt(2)  # cyl needs to be larger by sqrt(2) to contain cube
                obj['hax'] /= 2  # (/2 changes from full width to half)
            elif obj['type'].startswith('Sph'):
                obj['type'] = 1
            elif obj['type'].startswith('Ell'):
                obj['type'] = 1
            else:
                print('Unrecognized obj type')
                raise ValueError('Unrecognized obj type')

            obj['material'] = int(lns[18 + offset].strip().split()[1])
            obj['density'] = float(lns[19 + offset].strip().split()[1])
            offset += 21

        return obj


def euler_angs(R):
    p2 = np.real(np.arccos(R[2, 2]))
    if abs(np.sin(p2)) > 1e-6:
        p1 = np.arctan2(R[0, 2] / np.sin(p2), R[1, 2] / np.sin(p2))
        p3 = np.arctan2(R[2, 0] / np.sin(p2), -R[2, 1] / np.sin(p2))
    elif np.sign(np.cos(p2)) > 0:
        p2 = 0
        p3 = 0
        p1 = np.arctan2(R[0, 1], R[0, 0])
    else:
        p2 = np.pi
        p3 = 0
        p1 = np.arctan2(-R[0, 1], R[0, 0])

    Rnew = np.array([
        [np.cos(p1) * np.cos(p3) - np.sin(p1) * np.cos(p2) * np.sin(p3),
         np.cos(p1) * np.sin(p3) + np.sin(p1) * np.cos(p2) * np.cos(p3), np.sin(p1) * np.sin(p2)],
        [-np.sin(p1) * np.cos(p3) - np.cos(p1) * np.cos(p2) * np.sin(p3),
         -np.sin(p1) * np.sin(p3) + np.cos(p1) * np.cos(p2) * np.cos(p3), np.cos(p1) * np.sin(p2)],
        [np.sin(p2) * np.sin(p3), -np.sin(p2) * np.cos(p3), np.cos(p2)]
    ])

    err = np.linalg.norm(Rnew - R)
    if err > 1e-5:
        raise ValueError('Error computing Euler angles')

    return np.array([p1, p2, p3])


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
            line = file_contents[line_indices[i]+1 : line_indices[i+1]].decode('utf-8')
            lines.append(line)

    return lines