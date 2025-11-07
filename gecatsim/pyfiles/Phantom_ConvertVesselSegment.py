# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np
from gecatsim.pyfiles.euler_angles import euler_angles
from gecatsim.pyfiles.find_R import find_R

# Initialize parameters
std = np.array([1, 1, 1, 1, 1])
start_point_x = std * 0
start_point_y = np.arange(-1, 1.5, 0.5) * 100
start_point_z = std * (-65)

start_direction_t = std
start_direction_p = np.linspace(-0.4, 0.4, 5)
start_radius = std * 8

R_scaling_factor = np.concatenate((std * 0.7, std))
leng = np.concatenate((std * 3 * np.pi / 4, std * np.pi / 2))
ellipticity = np.concatenate((np.linspace(0.1, 0.1, 5), std * 0))
curvature = np.concatenate((std, std * 50))
torsion = np.concatenate((std * np.pi / 9, np.linspace(-np.pi / 3, -np.pi / 3, 5)))
tapering_nature = np.concatenate((std * 0, std * 0))

branch = np.array([0, 0, 0, 0, 0, 1, 2, 3, 4, 5])

n_objs = 0
endface = {}

objs = {
    'params': [],
    'str': [],
    'clip': []
}

for i in range(len(std) * 2):
    theta = leng[i]
    T = tapering_nature[i]
    delta = ellipticity[i]
    R1oR2 = R_scaling_factor[i]
    if branch[i]:
        R1 = endface[branch[i]]['R']
    else:
        R1 = start_radius[i]
    if delta == 0:
        C = curvature[i]
    else:
        C = 1

    phi = np.pi / 4 + T
    delta_1 = delta * np.cos(2 * (phi - theta / 2))
    delta_2 = delta * np.cos(2 * (phi + theta / 2))

    D_true = 1 / C ** 2 * (1 + delta / 2)
    E_true = 1 / C ** 2 * (1 - delta / 2)

    Da = 1 / C ** 2 * (1 + delta_1 / 2)
    Ea = 1 / C ** 2 * (1 + delta_2 / 2)

    alpha = 2 * Da
    beta = 2 * Ea

    if ((beta > alpha) and (R1oR2 > 1)) or ((beta < alpha) and (R1oR2 < 1)):
        phi = phi + np.pi / 2
        T = T + np.pi / 2
        delta_1 = delta * np.cos(2 * (phi - theta / 2))
        delta_2 = delta * np.cos(2 * (phi + theta / 2))

        D_true = 1 / C ** 2 * (1 + delta / 2)
        E_true = 1 / C ** 2 * (1 - delta / 2)

        Da = 1 / C ** 2 * (1 + delta_1 / 2)
        Ea = 1 / C ** 2 * (1 + delta_2 / 2)
        alpha = 2 * Da
        beta = 2 * Ea

    closest = np.sqrt(alpha / beta)
    if abs(1 - closest) > abs(1 - R1oR2):
        raise ValueError(f'critical R1oR2 is {closest:.5f}...decrease ellipticity, make R1oR2 farther from 1, or make T close to pi/4 (mod pi/2)')

    k = min([alpha, beta]) / 4
    calc = R1oR2 + 1
    step = k / 2
    c = 1
    if R1oR2 < 1:
        c = -c
    while abs(R1oR2 - calc) > 1e-10:
        s1 = np.sqrt(alpha - k + np.sqrt(alpha ** 2 - 2 * alpha * k))
        s2 = np.sqrt(alpha - k - np.sqrt(alpha ** 2 - 2 * alpha * k))
        s3 = np.sqrt(beta - k + np.sqrt(beta ** 2 - 2 * beta * k))
        s4 = np.sqrt(beta - k - np.sqrt(beta ** 2 - 2 * beta * k))
        calc = (s1 - s2) / (s3 - s4)
        err = R1oR2 - calc
        k = k + c * step * (2 * (err > 0) - 1)
        step = step / 2

    A = (s1 - s2) / 2 / R1
    if delta == 0:
        c = 0.001
        A = D_true
        calc = 2 * R1 * A + 1
        while abs(calc - 2 * R1 * A) > 1e-10:
            s1 = np.sqrt(alpha - k + np.sqrt(alpha ** 2 - 2 * alpha * k))
            s2 = np.sqrt(alpha - k - np.sqrt(alpha ** 2 - 2 * alpha * k))
            calc = (s1 - s2)
            k = max(0, min(min([alpha, beta]) / 2, k - c * (2 * R1 * A - calc)))
    s = k / A

    Sl = np.zeros((3, 3))
    Sr = np.zeros((3, 3))
    Sl[np.diag_indices(3)] = [A, A, A]
    Sr[0, 0] = D_true
    Sr[1, 1] = E_true

    delta_eff = delta * np.cos(2 * (phi - theta / 2))
    delta_eff_out = delta * np.cos(2 * (phi + theta / 2))
    a0t = np.array([np.cos(phi - theta / 2), np.sin(phi - theta / 2), 0]) * 1 / C / A * np.sqrt(1 + delta_eff / 2)
    a1t = np.array([np.cos(phi + theta / 2), np.sin(phi + theta / 2), 0]) * 1 / C / A * np.sqrt(1 + delta_eff_out / 2)
    cent_line_vec = np.array([-np.sin(phi - theta / 2), np.cos(phi - theta / 2), 0])
    cent_line_vec_out = np.array([-np.sin(phi + theta / 2), np.cos(phi + theta / 2), 0])
    curv_vec = np.array([-np.cos(phi - theta / 2), -np.sin(phi - theta / 2), 0])
    curv_vec_out = np.array([-np.cos(phi + theta / 2), -np.sin(phi + theta / 2), 0])
    if branch[i]:
        Dcent_line_vec = endface[branch[i]]['cent_line_vec']
        tvec2 = endface[branch[i]]['curv_vec']
        tvec1 = np.cross(tvec2, Dcent_line_vec)
        Dcurv_vec = np.cos(torsion[i]) * tvec1 + np.sin(torsion[i]) * tvec2
        start = endface[branch[i]]['center']
    else:
        Dcent_line_vec = np.array([np.sin(start_direction_t[i]) * np.cos(start_direction_p[i]),
                                   np.sin(start_direction_t[i]) * np.sin(start_direction_p[i]),
                                   np.cos(start_direction_t[i])])
        tvec2 = np.array([np.sin(start_direction_p[i]), -np.cos(start_direction_p[i]), 0])
        tvec1 = np.cross(tvec2, Dcent_line_vec)
        Dcurv_vec = np.cos(torsion[i]) * tvec1 + np.sin(torsion[i]) * tvec2
        start = np.array([start_point_x[i], start_point_y[i], start_point_z[i]])

    R = find_R(curv_vec, cent_line_vec, Dcurv_vec, Dcent_line_vec).T
    center = start - R.T @ a0t
    Dcent_line_vec_out = R.T @ cent_line_vec_out
    Dcurv_vec_out = R.T @ curv_vec_out

    clip = np.column_stack((-Dcent_line_vec, Dcent_line_vec_out))
    clip = np.vstack((clip, (center.T @ clip).reshape(1, -1)))

    euler_angs = euler_angles(R) * 180 / np.pi

    endface[i] = {
        'center': center + R @ a1t,
        'R': R1 / R1oR2,
        'curv_vec': Dcurv_vec_out,
        'cent_line_vec': Dcent_line_vec_out
    }

    j = i + n_objs
    inten = 1
    trans = 0
    objs['params'].append(np.concatenate((center, [k, delta, A], euler_angs, [inten, 7, trans, C, theta, T])))
    objs['str'].append('')
    objs['clip'].append(clip)

print("Clip for segment 1:")
print(objs['clip'][0])

print("Clip for segment 6:")
print(objs['clip'][5])

