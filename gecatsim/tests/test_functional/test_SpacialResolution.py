import unittest.mock

import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft2, fftshift
from scipy.interpolate import interpn

import gecatsim as xc
import gecatsim.reconstruction.pyfiles.recon as recon


class Test_Functional_SpacialResolution(unittest.TestCase):

    def test_create_spatialresolution(self):

        ##--------- Initialize
        ct = xc.CatSim("../examples/mtf_examples/Phantom_Sample_Analytic", "../examples/mtf_examples/Scanner_Sample_generic",
                       "../examples/mtf_examples/Protocol_Sample_axial")  # initialization

        ##--------- Make changes to parameters (optional)
        ct.resultsName = "W_wire"

        ct.phantom.filename = 'tungsten_wire.ppm'
        ct.phantom.centerOffset = [30, 50, 0]

        ct.protocol.viewsPerRotation = 1000
        ct.protocol.viewCount = ct.protocol.viewsPerRotation
        ct.protocol.stopViewId = ct.protocol.viewCount - 1

        ct.scanner.detectorRowsPerMod = 1
        ct.scanner.detectorRowCount = ct.scanner.detectorRowsPerMod

        ct.physics.colSampleCount = 12
        ct.physics.rowSampleCount = 2
        ct.physics.srcXSampleCount = 12
        ct.physics.srcYSampleCount = 2
        ct.physics.viewSampleCount = 6

        ##--------- Run simulation
        ct.run_all()  # run the scans defined by protocol.scanTypes

        ##--------- Reconstruction
        ct.recon.fov = 50.0
        ct.recon.sliceCount = ct.scanner.detectorRowCount  # number of slices to reconstruct
        ct.recon.sliceThickness = 0.568  # reconstruction inter-slice interval (in mm)
        ct.recon.centerOffset = -np.array(ct.phantom.centerOffset)
        ct.do_Recon = 1
        recon.recon(ct)

        imgFname = "%s_%dx%dx%d.raw" % (ct.resultsName, ct.recon.imageSize, ct.recon.imageSize, ct.recon.sliceCount)
        fig_name = imgFname + '.png'
        img = xc.rawread(imgFname, [ct.recon.sliceCount, ct.recon.imageSize, ct.recon.imageSize], 'float')
        plt.imshow(img[np.int32(ct.recon.sliceCount / 2), :, :], cmap='gray', vmin=-1000, vmax=200)
        plt.savefig(fig_name)

        assert self.getFileSize(
            "%s.air" % (ct.resultsName)) == ct.scanner.detectorRowCount * ct.scanner.detectorColCount
        assert self.getFileSize(
            "%s.offset" % (ct.resultsName)) == ct.scanner.detectorRowCount * ct.scanner.detectorColCount
        assert self.getFileSize(
            "%s.prep" % (
                ct.resultsName)) == ct.scanner.detectorRowCount * ct.scanner.detectorColCount * ct.protocol.viewsPerRotation
        assert self.getFileSize(
            imgFname) == ct.recon.imageSize * ct.recon.imageSize * ct.recon.sliceCount

    def test_functional_spatialresolution(self):
        # Main script
        np_val = 512
        img = self.rawread('./W_wire_512x512x1.raw', (np_val, np_val), 'float32')

        fov = 50  # mm
        cx = np_val / 2 + 0.5
        cy = np_val / 2 + 0.5

        # Detrend and normalize
        roi = ~self.fovimg(np_val, np_val, 1, 50, cx, cy)
        meanBkg = np.mean(img[roi])
        img = img - meanBkg
        img = img / np.max(img)

        # Select ROI
        threshold = 0.1
        tmp = img.copy()
        tmp[tmp >= threshold] = 1
        tmp[tmp < threshold] = 0
        dx = np.sum(np.sum(tmp, axis=0) > 0)
        dy = np.sum(np.sum(tmp, axis=1) > 0)
        r_roi = max(dx, dy) / 2
        roi = self.fovimg(np_val, np_val, 1, r_roi, cx, cy)
        img = img * roi

        # 2D FFT
        otf = fft2(img)
        otf = np.abs(otf) / np.abs(otf[0, 0])
        otf = fftshift(otf)
        otf = np.real(otf)

        # Radial mean
        print('Started calculating MTF...')

        radial_profile, xx = self.RadialProfile(otf)
        faxis = xx / fov * 10  # lp/cm
        cent = len(faxis) // 2
        faxis = faxis[cent:]
        print('Still calculating MTF...')
        radial_profile = (radial_profile[cent::-1] + radial_profile[cent:]) / 2
        mtf_curve = radial_profile / radial_profile[0]
        print('Still calculating MTF...Almost done...')
        mtf50 = faxis[np.where(mtf_curve >= 0.5)[0][-1]]
        mtf10 = faxis[np.where(mtf_curve >= 0.1)[0][-1]]
        mtf5 = faxis[np.where(mtf_curve >= 0.05)[0][-1]]

        print(f'MTF (0.50, 0.10, 0.05) = {mtf50:.1f}, {mtf10:.1f}, {mtf5:.1f}')

        assert f'{mtf50:.1f}' == '4.6'
        assert f'{mtf10:.1f}'== '7.8'
        assert f'{mtf5:.1f}' == '8.6'


    def getFileSize(self, fileName):
        fileDataRaw = np.fromfile(fileName, "<f")
        return int(len(fileDataRaw))

    def rawread(self, filename, shape, dtype):
        with open(filename, 'rb') as f:
            data = np.fromfile(f, dtype=dtype)
        return data.reshape(shape)

    def fovimg(self, np1, np2, val1, val2, cx, cy):
        x = np.arange(np1) - cx
        y = np.arange(np2) - cy
        X, Y = np.meshgrid(x, y)
        return np.sqrt(X ** 2 + Y ** 2) <= val2

    def RadialProfile(self, img):
        nAngSample = 100
        n = img.shape[0]
        v = np.arange(1, n + 1) - (n / 2 + 0.5)
        X, Y = np.meshgrid(v, v)

        if n % 2 == 1:
            xx = v
        else:
            xx = np.arange(-(n // 2 - 1), n // 2)

        radial_profile = np.zeros_like(xx, dtype=float)
        for ii in range(nAngSample):
            ang = (ii * np.pi) / nAngSample
            x1 = xx * np.cos(ang)
            y1 = xx * np.sin(ang)
            points = (v, v)
            xi = np.vstack((x1, y1)).T
            p = interpn(points, img, xi, method='cubic')
            radial_profile += p / nAngSample

        return radial_profile, xx
