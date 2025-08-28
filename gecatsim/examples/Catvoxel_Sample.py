import gecatsim as xc
from gecatsim.pyfiles import catvoxel
'''
examples to use catvoxel.
Jiayong Zhang
'''

ct = xc.CatSim()
cfg = ct.get_current_cfg()
# -------------------
# some common settings

# Make separate volumes of densities (in gm/cm^3) for each materia
cfg.material_volumes = 1
# oversampling in making volumes
cfg.vol_os = 11
# Sets the kV when cfg.material_volumes is set to 0
# in this case, the attenuation coefficients at this kVp will be saved
# cfg.make_img_kv = 120
# save results
cfg.write_vp = 1

cfg.recon.imageSize = 512
cfg.recon.sliceCount = 10
cfg.recon.sliceThickness = 0.568
cfg.phantom.scale = 1.0

# analytic phantom
print("running catvoxel for analytic phantom...")
cfg.phantom.projectorCallback = 'C_Projector_Analytic'
cfg.phantom.callback = 'Phantom_Analytic'
cfg.phantom.filename = 'CTDI_16cm_WaterAirPEBoneChambers.ppm'
cfg.recon.fov = 500.0
cfg.recon.centerOffset = [0, 0, 0]
catvoxel(cfg)

# NCAT phantom
print("running catvoxel for NCAT phantom...")
cfg.phantom.projectorCallback = 'C_Projector_NCAT'
cfg.phantom.callback = 'Phantom_NCAT'
cfg.phantom.filename = 'vmale50_chest_less_surfaces.nrb'
cfg.recon.fov = 400.0
cfg.recon.centerOffset = [0, 0, 0]
catvoxel(cfg)

# polygonal phantom
print("running catvoxel for polygonal phantom...")
cfg.phantom.projectorCallback = 'C_Projector_Polygon'
cfg.phantom.callback = 'Phantom_Polygonal'
cfg.phantom.filename = 'female_adult_average_lung_lesions_reduced.nrb'
cfg.recon.fov = 50.0
cfg.recon.centerOffset = [-57, -51, 102]
catvoxel(cfg)
