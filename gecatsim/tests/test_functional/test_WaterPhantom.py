import unittest.mock
from unittest.mock import patch

import os
import numpy as np

import gecatsim as xc
import gecatsim.reconstruction.pyfiles.recon as recon
from gecatsim.pyfiles.CommonTools import *
#from gecatsim.pyfiles.Detection_PC import Detection_PC

import matplotlib.pyplot as plt
from operator import __getitem__

class Test_Functional_WaterPhantom(unittest.TestCase):

    def test_create_waterphantom(self):

        #ct = xc.CatSim( "all_in_one_config_generic_scanner") #initalization
        ct = xc.CatSim("../examples/vct_examples/Phantom_Sample_Analytic", 
               "../examples/vct_examples/Physics_Sample", 
              "../examples/vct_examples/Protocol_Sample_axial", 
              "../examples/vct_examples/Recon_Sample_2d", 
              "../examples/vct_examples/Scanner_Sample_generic") #initalization
         
        # Make changes to parameters as required
        ct.resultsName = "Water_IQ_test"
        ct.phantom.filename ='W20.ppm'  # water phantom filename
        ct.recon.fov = 300.0           # diameter of the reconstruction field-of-view (in mm)
        ct.recon.sliceCount = 4        # number of slices to reconstruct
        ct.physics.energyCount = 24
        ct.physics.monochromatic = -1
        ct.physics.colSampleCount = 2
        ct.physics.rowSampleCount = 2
        ct.physics.srcXSampleCount = 2
        ct.physics.srcYSampleCount = 2
        ct.physics.viewSampleCount = 1
        ct.protocol.mA = 130                           # tube current (in mA)
        ct.protocol.spectrumFilename = "xcist_kVp120_tar7_bin1.dat" # name of the spectrum file

        ct.physics.callback_post_log = 'Prep_BHC_Accurate'
        ct.physics.EffectiveMu = 0.2
        ct.physics.BHC_poly_order = 5
        ct.physics.BHC_max_length_mm = 300
        ct.physics.BHC_length_step_mm = 10

        # Run simulation
        ct.run_all()
         
        #assert(os.path.exists('all_in_one_config_generic_scanner.cfg') == True)
        #assert(os.path.exists('all_in_one_config.cfg') == True)

        # Reconstruction
        ct.do_Recon = 1
        recon.recon(ct)
         
        # Show results
        imgFname = "%s_%dx%dx%d.raw" %(ct.resultsName, ct.recon.imageSize, ct.recon.imageSize, ct.recon.sliceCount)
        img = xc.rawread(imgFname, [ct.recon.sliceCount, ct.recon.imageSize, ct.recon.imageSize], 'float')
        plt.imshow(img[2,:,:], cmap='gray', vmin=-200, vmax=200)
        plt.savefig(imgFname)
        #plt.show()

    def test_functional_waterphantom(self):

        fileName = "./Water_IQ_test_512x512x4.raw"
        imgRows = 512
        imgCols = 512
        imgDataRaw = np.fromfile(fileName, "<f")
        numImgs = int(len(imgDataRaw) / (imgRows * imgCols))
        imgStackRaw = np.reshape(imgDataRaw, (numImgs, imgRows, imgCols))
        imgNum = 1
        plt.imshow(imgStackRaw[imgNum, :, :], cmap='gray')
               
        #Select ROI at  center of the phantom and get HU value and noise value

        roiMatrix_C = imgStackRaw[imgNum, 236:276,236 :276]
        meanHU_C = np.mean(roiMatrix_C)
        noise_ROI_40_40=np.std(roiMatrix_C)
        CT_Number_Water= 0
        CT_num_Deviation=CT_Number_Water - meanHU_C

        print("CT_Number_Water at center:", round(meanHU_C,2),"HU")

        if meanHU_C <3 and meanHU_C>-3:
            print('CT number is wthin the range, test pass!')
        else:
            print('CT number is not within the range, test failed!')

        print("Noise_at centre:" ,round(noise_ROI_40_40,2))
         
        
            
        roiMatrix_P1 = imgStackRaw[imgNum, 105:145,245 :285]
        roiMatrix_P2 = imgStackRaw[imgNum, 365:405,245 :285]
        roiMatrix_P3 = imgStackRaw[imgNum, 245:285,105 :145]
        roiMatrix_P4 = imgStackRaw[imgNum, 245:285,365 :405]
           
        #calculate mean deviation from centre value
        meanHU_P1 = np.mean(roiMatrix_P1)
        meanHU_P2 = np.mean(roiMatrix_P2)
        meanHU_P3 = np.mean(roiMatrix_P3)
        meanHU_P4 = np.mean(roiMatrix_P4)

        deviation= [meanHU_P1-meanHU_C,meanHU_P2-meanHU_C,meanHU_P3-meanHU_C,meanHU_P4-meanHU_C]
        uniformity=sum(deviation)/len(deviation)

        print('CT number at periphery: ',round(meanHU_P1,2),',',round(meanHU_P2,2),',',round(meanHU_P3,2),',',round(meanHU_P4,2))
        if meanHU_P1 <3 and meanHU_P1 >-3:
            if meanHU_P2 <3 and meanHU_P2 >-3:
                if meanHU_P3 <3 and meanHU_P3 >-3:
                    if meanHU_P4 <3 and meanHU_P4 >-3:
                        print( 'CT number at periphery is within the range, test pass!')
        else: 
            print ('CT number at periphery is not within the range, test failed!')

        print("CT Number Uniformity:",round(uniformity,2),"HU")

        #calsulate Noise uniformity
        noise_P1 = np.std(roiMatrix_P1)
        noise_P2 = np.std(roiMatrix_P2)
        noise_P3 = np.std(roiMatrix_P3)
        noise_P4 = np.std(roiMatrix_P4)

        Noise_deviation= [noise_P1-noise_ROI_40_40,noise_P2-noise_ROI_40_40,noise_P3-noise_ROI_40_40,noise_P4-noise_ROI_40_40]
        Noise_uniformity=sum(Noise_deviation)/len(Noise_deviation)

        print('Noise at periphery:', round(noise_P1,2),',',round(noise_P2,2),',',round(noise_P3,2),',',round(noise_P4,2))
        print("Noise Uniformity:",round(Noise_uniformity,2))

        #plt.show()

if __name__ == "__main__":
    unittest.main()
