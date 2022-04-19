import os
import copy
import numpy as np
from catsim.CommonTools import rawread
from catsim.CatSim import CatSim 

'''
prep: viewxrowxcol
'''
# let's first try with a single view and compare it with matlab
def Prep_BHC_Accurate(cfg, view_in, view_idx):
    ct = CatSim()
    ct.cfg_to_self(cfg)
    # if BHC_vec_fname already exists
    if view_idx == ct.protocol.startViewId:
        if hasattr(ct.physics, "BHC_vec_fname"):
            if os.path.isfile(ct.physics.BHC_vec_fname):
                poly_coef = np.load(ct.physic.BHC_vec_fname)
            else:
                poly_coef = gen_BHC_vec(ct)
                np.save(ct.physics.BHC_vec_fname, poly_coef)
        else:
            poly_coef = gen_BHC_vec(ct)

    breakpoint()
    view_out = poly_coef[:, 0]
    num_poly_coefs = ct.physics.BHC_poly_order+1
    for i in range(1, num_poly_coefs):
        view_out = view_out*view_in + poly_coef[:,i]

    #breakpoint()
    return view_out

def gen_BHC_vec(ct):

    poly_order = ct.physics.BHC_poly_order
    num_poly_coefs = poly_order+1

    if not hasattr(ct.physics, "BHC_material") or ct.physics.BHC_material == "":
        mt = "water"
    else:
        mt = ct.physics.BHC_material

    max_length_mm = ct.physics.BHC_max_length_mm
    length_step_mm = ct.physics.BHC_length_step_mm
    num_length = int(max_length_mm/length_step_mm)
    #preadjust = []#TODO: how to do preadjust in python; I dont' think that we need to worry about this now

    # does not work now, will test in a full run
    #if hasattr(cfg.det, "totalNumCells"):
    #    ndet = cfg.det.totalNumCells
    #else:
    #    ndet = cfg.scanner.detectorColCount*cfg.scanner.detectorRowCount
    ndet = ct.scanner.detectorColCount*ct.scanner.detectorRowCount
    sig = np.zeros((num_length, ndet))

    # perform air scan
    airscan_ct = copy.deepcopy(ct)
    #airscan_cfg.sim.thisScanType = [1, 0, 0]
    airscan_ct.protocol.scanTypes = [1, 0, 0]
    airscan_ct.protocol.airViewCount = 1
    #airscan_ct.protocol.viewCount = 1
    #breakpoint()
    airscan_ct.physics.recalcFilt = 1
    #airscan_ct.protocol.bowtie = []
    #airscan_ct.protocol.flatFilter = []
    #airscan_ct.scanner.detectorPrefilter = []
    # TODO: need to find what is in cxist for varialbre_kv_ma, and recompute_prefilter
    airscan_ct.air_scan()
    #air_scan_cfg = airscan_ct.self_to_cfg()
    mtViews = rawread(airscan_ct.resultsName+".air", [ndet, 1], 'float')
    mtViews = mtViews.reshape(airscan_ct.scanner.detectorColCount, airscan_ct.scanner.detectorRowCount)
    I0 = mtViews.T
    #breakpoint()
    #assert np.allclose(I0, mat_I0, atol=1.E-6), 'Error! I02 do not match'
    #breakpoint()

    # now calculate 
    # need to start from 1 otherwise thick will be 0
    for j in range(num_length):
        thick = (1+j)*length_step_mm

        airscan_ct = copy.deepcopy(ct)
        airscan_ct.protocol.scanTypes = [1, 0, 0]
        airscan_ct.protocol.viewCount = 1
        #airscan_ct.protocol.bowtie = []
        #airscan_ct.protocol.flatFilter = []
        #airscan_ct.scanner.detectorPrefilter = []
        airscan_ct.scanner.detectorPrefilter += [mt, thick]
        airscan_ct.air_scan()
        rawViews = rawread(airscan_ct.resultsName+".air", [ndet, 1], 'float')
        rawViews = rawViews.reshape(airscan_ct.scanner.detectorColCount, airscan_ct.scanner.detectorRowCount)
        sig[j] = -np.log(rawViews.T/I0)

    # final run to save .scan
    airscan_ct = copy.deepcopy(ct)
    airscan_ct.protocol.scanTypes = [1, 0, 0]
    #airscan_ct.protocol.bowtie = []
    #airscan_ct.protocol.flatFilter = []
    #airscan_ct.scanner.detectorPrefilter = []
    airscan_ct.air_scan()

    # correct and fit signal
    dessig = cfg.physics.EffectiveMu*length_step_mm/10*(1+np.arange(num_length))
    poly_coef = np.zeros((ndet, num_poly_coefs))
    for n in range(ndet):
        poly_coef[n] = np.polyfit(sig[:,n], dessig, poly_order)

    breakpoint()
    return poly_coef

if __name__=="__main__":
    # currently there is no protocol.maxPrep in this example, thus we can use prep
    # we may need to use .scan, .air, .offset to benchmark in some cases
    # actually maxPrep is 9 in calculations
    import catsim as xc
    #prep = xc.rawread("FB_head_from_matlab.prep", [10,1,888], 'float')
    ct = xc.CatSim("Phantom_Sample","Physics_Sample", "Scanner_Sample_generic", "Protocol_Sample_axial")
    #breakpoint()
    ct.resultsName = "my_test"
    cfg = ct.self_to_cfg()

    from scipy import io as sio
    mat_in = sio.loadmat("Prep_BHC_Accurate_view_in.mat")['view_in']
    mat_in = mat_in[:,0]
    mat_out = sio.loadmat("Prep_BHC_Accurate_view_out.mat")['view_out']
    mat_I0 = sio.loadmat("Prep_BHC_Accurate_I0.mat")['I0']
    corr_prep = Prep_BHC_Accurate(cfg, mat_in, 0)
    try:
        assert np.allclose(corr_prep, mat_out, atol=1.E-8), 'Error! corr_prep does not match mat_out'
    except AssertionError as err:
        print(err)
        breakpoint()
