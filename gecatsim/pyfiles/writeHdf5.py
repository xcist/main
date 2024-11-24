# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import os
import numpy as np
import h5py
from gecatsim.pyfiles.readViewWeighting import readViewWeighting

def writeHdf5(cfg, ofolder):
    row_thickness_mm = cfg.recon_slice_thickness
    if not ofolder:
        ofolder = "."
    if not os.path.exists(ofolder):
        raise ValueError("Folder does not exist in the path: specify a foldername using the second input argument")

    print(f"Creating an ADF at location {ofolder}")
    scalefactor = 1000.0 / (cfg.recon_mu * 0.1) / (4.0 * cfg.col_size * cfg.sid / cfg.sdd)
    infile = f"{cfg.results_basename}.prep"
    fp = open(infile, "rb")
    total_n_views = cfg.total_n_viewsverbo

    ofile = f"{cfg.results_basename}.h5"
    if os.path.isfile(ofile):
        os.remove(ofile)

    with h5py.File(ofile, "w") as h5f:
        for i in range(1, total_n_views + 1):
            sino = np.fromfile(fp, dtype=np.float32, count=cfg.col_count * cfg.row_count)
            sino = sino.reshape((cfg.col_count, cfg.row_count))
            sino *= scalefactor
            h5f.create_dataset(f"/counts/{i}", data=sino)

        # Add other datasets and attributes...
        h5f.create_dataset("/info/gantry_angle_rad", data=cfg.betas)
        h5f.create_dataset("/info/current_ma", data=np.repeat(cfg.mA, total_n_views))
        h5f.create_dataset("/info/integration_period_s",
                           data=np.repeat(cfg.rotation_period / cfg.views_per_rotation, total_n_views))
        h5f.create_dataset("/info/voltage_kv", data=np.repeat(120, total_n_views))

        # Add attributes
        h5f.attrs["det_angle_rad"] = (np.pi / 180) * (cfg.gammas[0, 1] - cfg.gammas[0, 0])
        h5f.attrs["rot_time_sec"] = cfg.rotation_period
        h5f.attrs["sbd"] = cfg.sbd
        h5f.attrs["sdd"] = cfg.sdd
        h5f.attrs["sid"] = cfg.sid
        center_detector = (cfg.col_count + 1) / 2 - cfg.col_offset - 1
        h5f.attrs["cellNumberAtTheta"] = center_detector
        trigger_freq_hz = cfg.views_per_rotation / cfg.rotation_period
        h5f.attrs["trigger_freq_hz"] = trigger_freq_hz
        h5f.attrs["views_per_rotation"] = cfg.views_per_rotation
        h5f.attrs["kvp"] = 120
        h5f.attrs["mamp"] = cfg.mA
        h5f.attrs["row_at_iso_mm"] = row_thickness_mm

        # Add scan related parameters
        imgTable = calcImgTableParams(cfg, row_thickness_mm)
        h5f.attrs["pitch_per_view_mm"] = imgTable.pitch_per_view_mm
        h5f.attrs["pitch_rows"] = imgTable.pitch_rows
        h5f.attrs["scan_pitch_ratio"] = imgTable.scan_pitch_ratio
        h5f.attrs["first_image_location_mm"] = imgTable.first_image_location_mm
        h5f.attrs["number_of_images"] = imgTable.number_of_images
        h5f.attrs["first_view_for_image"] = imgTable.first_view_for_image
        h5f.attrs["number_of_views_for_image"] = imgTable.number_of_views_for_image
        h5f.attrs["image_location_increment_mm"] = imgTable.image_location_increment_mm
        h5f.attrs["table_speed_mm_sec"] = imgTable.table_speed_mm_sec
        h5f.attrs["actual_table_end_mm"] = imgTable.actual_table_end_mm
        h5f.attrs["actual_table_start_mm"] = imgTable.actual_table_start_mm
        h5f.attrs["table_direction"] = imgTable.table_direction
        h5f.attrs["number_of_views_per_image"] = imgTable.number_of_views_per_image

        # Add view weighting parameters for the image table module
        if cfg.table_speed > 0:
            collimation = cfg.row_count
            pitch = imgTable.pitch_rows
            vparams = readViewWeighting(collimation, pitch)
        else:
            vparams = {
                "beta_0": -1.0,
                "beta_t": -1.0,
                "vct_k": 0.0,
                "vct_q": 0.0,
                "vct_r": 0.0,
                "zsf": 1.0,
                "kw": 0.0
            }

        h5f.attrs["vw_vct_beta_0"] = vparams["beta_0"]
        h5f.attrs["vw_vct_beta_t"] = vparams["beta_t"]
        h5f.attrs["vw_vct_k"] = vparams["vct_k"]
        h5f.attrs["vw_vct_q"] = vparams["vct_q"]
        h5f.attrs["vw_vct_r"] = vparams["vct_r"]
        h5f.attrs["thick_factor"] = vparams["zsf"]
        h5f.attrs["k_weight"] = vparams["kw"]

        # MW: SKIP calling Recon_FilterParams_Get()
        # Add filter module in the info header
        # fparams = Recon_FilterParams_Get(cfg.recon_filter, True)
        # h5f.attrs["filter_typ"] = cfg.recon_filter
        # h5f.attrs["window_type"] = np.int16(fparams.window)
        # h5f.attrs["freq_interp"] = np.int16(fparams.freqInterp)
        # h5f.attrs["fmax"] = fparams.fmax
        # h5f.attrs["fpixf"] = fparams.fpix
        # h5f.attrs["flow"] = fparams.flow
        # h5f.attrs["poly_coeff"] = fparams.p

    print(f"Saved HDF5 file: {ofile}")


def getNumExtraViews(cfg):
    centerDet = (1.0 + cfg.col_count) * 0.5 - cfg.col_offset - 1.0
    gamma = cfg.col_size * centerDet / cfg.sdd
    dBeta = 2.0 * np.pi / cfg.views_per_rotation
    LAGR_BETA_INTERP_ORDER = 6
    nExtraViews = int(np.floor(gamma / dBeta) + LAGR_BETA_INTERP_ORDER / 2)
    return nExtraViews


def calcImgTableParams(cfg, row_thickness_mm):
    imgTable = {}

    if cfg.table_speed > 0:
        pitch_rows = round(cfg.table_speed * cfg.rotation_period / row_thickness_mm)
        relPitch = pitch_rows / cfg.row_count
        normpitch = abs(relPitch)

        if normpitch < 0.7:
            nViewsForImage = 2 * cfg.views_per_rotation
        else:
            nViewsForImage = cfg.views_per_rotation

        nExtraViews = getNumExtraViews(cfg)
        firstCenterView = cfg.start_view + nExtraViews + nViewsForImage // 2
        incrPerView = cfg.table_speed * cfg.rotation_period / cfg.views_per_rotation
        firstPossibleImgPos = cfg.start_z + incrPerView * firstCenterView
        lastCenterView = cfg.start_view + cfg.total_n_views - 1 - nExtraViews - nViewsForImage // 2
        lastPossibleImgPos = cfg.start_z + incrPerView * lastCenterView

        firstPossibleReqImage = np.ceil((firstPossibleImgPos - cfg.recon_zcenter) / cfg.recon_slice_thickness)
        lastPossibleReqImage = np.floor((lastPossibleImgPos - cfg.recon_zcenter) / cfg.recon_slice_thickness)

        nImages = 0
        imgPos = cfg.recon_zcenter

        if firstPossibleReqImage > lastPossibleReqImage:
            print("WARNING: helicalImageTable could not be built for first image at", cfg.recon_zcenter,
                  "mm and last image at", cfg.recon_zcenter + (cfg.recon_planes - 1) * cfg.recon_slice_thickness, "mm")
        else:
            nImages = int(lastPossibleReqImage - firstPossibleReqImage + 1)
            imgPos = cfg.recon_zcenter + firstPossibleReqImage * cfg.recon_slice_thickness

            if firstPossibleReqImage > 0:
                print("WARNING: helicalImageTable: first image position will be at", imgPos,
                      "instead of", cfg.recon_zcenter)
            if lastPossibleReqImage < cfg.recon_planes - 1:
                print("WARNING: helicalImageTable: last image position will be at", imgPos,
                      "instead of", cfg.recon_zcenter + lastPossibleReqImage * cfg.recon_slice_thickness,
                      "mm")

        print("helicalImageTable: number of images is", nImages, "(requested", cfg.recon_planes, ")")
        print("First possble image position is", imgPos, "mm, last is",
              imgPos + (nImages - 1) * cfg.recon_slice_thickness, "mm")

    else:
        incrPerView = 0
        pitch_rows = 0
        normpitch = 0.0
        imgPos = 0 - cfg.recon_slice_thickness * (cfg.recon_planes // 2)
        nImages = cfg.recon_planes
        nViewsForImage = cfg.total_n_views

    imgTable["pitch_per_view_mm"] = incrPerView
    imgTable["pitch_rows"] = pitch_rows
    imgTable["scan_pitch_ratio"] = f"{normpitch:.3f}"
    imgTable["first_image_location_mm"] = imgPos
    imgTable["number_of_images"] = nImages
    imgTable["first_view_for_image"] = cfg.start_view
    imgTable["number_of_views_for_image"] = cfg.total_n_views
    imgTable["image_location_increment_mm"] = cfg.recon_slice_thickness
    imgTable["table_direction"] = np.sign(cfg.table_speed)
    imgTable["actual_table_start_mm"] = cfg.start_z
    imgTable["actual_table_end_mm"] = cfg.start_z + incrPerView * cfg.total_n_views
    imgTable["table_speed_mm_sec"] = cfg.table_speed
    imgTable["number_of_views_per_image"] = nViewsForImage

    return imgTable

