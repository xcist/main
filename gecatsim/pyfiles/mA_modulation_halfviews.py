# Copyright 2025, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

# -----------------------------------------------------------------------
#   Program Name: mA_modulation_halfviews.py
#   Author:  Samit Basu, Bruno De Man, Jed Pack (GE Global Research)
#   Organization: GE Global Research and GE Healthcare, General Electric Company
#   Version:  6.0.3
#   Date:  Feb 3, 2015
#   Class: GE Confidential. General Electric Proprietary Data (c) 2012 General Electric Company
#
# History:
#   2020-03-19 Debashish Pal (GE Healthcare)
#              New file to simulate half views getting higher mA than other
#              half views.
# -----------------------------------------------------------------------

import matplotlib.pyplot as plt

def mA_modulation_halfviews(cfg, viewnr):

    print(f'Calculating mA modulation using {__file__}')

    st = cfg.start_view + cfg.views_per_rotation // 4
    range_views = list(range(st, st + cfg.views_per_rotation // 2))

    if viewnr in range_views:
        mA_scale_factor = 2
    else:
        mA_scale_factor = 1

    if viewnr == 1:
        plt.figure(1)
        plt.scatter(viewnr, mA_scale_factor, color='r')
        plt.ion()  # Enable interactive mode

    plt.scatter(viewnr, mA_scale_factor, color='r')
    plt.show()

    return mA_scale_factor