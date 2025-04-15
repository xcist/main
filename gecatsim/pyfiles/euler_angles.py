# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import numpy as np

def euler_angles(R):
    print(f'Calculating Euler angles using {__file__}')

    p2 = np.real(np.arccos(R[2, 2]))  # could also be -acos
    if abs(np.sin(p2)) > 1e-6:
        p1 = np.arcsin(R[0, 2] / np.sin(p2))  # could also be pi-asin
        if np.sign(np.cos(p1) * np.sin(p2)) != np.sign(R[1, 2]):
            p1 = np.pi - p1
        p3 = np.real(np.arcsin(R[2, 0] / np.sin(p2)))  # could also be pi-asin
        if np.sign(-np.cos(p3) * np.sin(p2)) != np.sign(R[2, 1]):
            p3 = np.pi - p3
    elif np.sign(np.cos(p2)) > 0:
        p2 = 0
        p3 = 0
        p1 = np.arccos(R[0, 0])
        if np.sign(-np.sin(p1)) != np.sign(R[1, 0]):
            p1 = -p1
    else:
        p2 = np.pi
        p3 = 0
        p1 = np.arccos(R[0, 0])
        if np.sign(-np.sin(p1)) != np.sign(R[1, 0]):
            p1 = -p1

    Rnew = np.array([
        [np.cos(p1) * np.cos(p3) - np.sin(p1) * np.cos(p2) * np.sin(p3), np.cos(p1) * np.sin(p3) + np.sin(p1) * np.cos(p2) * np.cos(p3), np.sin(p1) * np.sin(p2)],
        [-np.sin(p1) * np.cos(p3) - np.cos(p1) * np.cos(p2) * np.sin(p3), -np.sin(p1) * np.sin(p3) + np.cos(p1) * np.cos(p2) * np.cos(p3), np.cos(p1) * np.sin(p2)],
        [np.sin(p2) * np.sin(p3), -np.sin(p2) * np.cos(p3), np.cos(p2)]
    ])
    err = np.linalg.norm(Rnew - R)
    if err > 1e-5:
        raise ValueError('Error computing Euler angles')
    ea = np.array([p1, p2, p3])
    return ea