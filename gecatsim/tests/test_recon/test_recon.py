from gecatsim.pyfiles.CommonTools import *
import gecatsim.reconstruction.pyfiles.recon as recon
import gecatsim as xc
import unittest.mock
from unittest.mock import patch
import sys

class Test_recon(unittest.TestCase):
    def test_recon_do_recon(self):

        ct = self.initialize()
        retval = recon.recon(ct)

        assert(sys.getsizeof(retval)>0)
        assert(ct==retval)

    @patch('gecatsim.reconstruction.pyfiles.recon.xc.rawwrite', create=True)
    def test_recon_do_recon_save_img_vol(self, recon_xc_raw_write):

        ct = self.initialize()
        ct.recon.saveImageVolume = True
        retval = recon.recon(ct)

        imageVolume3D_size_string = str(ct.recon.imageSize) + 'x' + str(ct.recon.imageSize) + 'x' + str(
            ct.recon.sliceCount)
        fname = ct.resultsName + '_' + imageVolume3D_size_string + '.raw'

        assert recon_xc_raw_write.call_count == 1

    def test_recon_load_recon(self):

        ct = self.initialize()
        ct.do_Recon = False
        retval = recon.recon(ct)

        assert(sys.getsizeof(retval)>0)
        assert(ct==retval)

    @patch('gecatsim.reconstruction.pyfiles.recon.xc.rawwrite', create=True)
    def test_recon_do_recon_units(self, recon_xc_raw_write):

        ct = self.initialize()

        prep = recon.load_prep(ct)
        # cfg.recon.reconType is the recon function's name
        imageVolume3D = feval("gecatsim.reconstruction.pyfiles." + ct.recon.reconType, ct, prep)

        imageVolume3D_size_string = str(ct.recon.imageSize) + 'x' + str(ct.recon.imageSize) + 'x' + str(
            ct.recon.sliceCount)
        fname = ct.resultsName + '_' + imageVolume3D_size_string + '.raw'
        imageVolume3D = imageVolume3D.transpose(2, 0, 1)
        imageVolume3D = imageVolume3D.copy(order='C')

        ct.recon.unit = '/mm'
        retval = recon.recon(ct)

        assert (recon_xc_raw_write.call_args.args[1] == imageVolume3D).all()

        imageVolume3D = imageVolume3D * 10

        ct.recon.unit = '/cm'
        retval = recon.recon(ct)

        assert (recon_xc_raw_write.call_args.args[1] == imageVolume3D).all()

    @patch('gecatsim.reconstruction.pyfiles.recon.xc.rawwrite', create=True)
    def test_recon_save_single_images(self, recon_xc_raw_write):

        ct = self.initialize()
        ct.do_Recon = False
        ct.recon.saveSingleImages = True
        retval = recon.recon(ct)

        assert recon_xc_raw_write.call_count == ct.recon.sliceCount

    @patch('gecatsim.reconstruction.pyfiles.recon.plt.figure', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.imshow', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.title', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.close', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.draw', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.pause', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.savefig', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.axis', create=True)
    def test_recon_image_pictures_display(self, recon_plt_axis, recon_plt_savefig, recon_plt_pause, recon_plt_draw, recon_plt_close, recon_plt_title, recon_plt_imshow, recon_plt_figure):
        ct = self.initialize()
        ct.do_Recon = False
        ct.recon.displayImagePictures = True
        ct.waitForKeypress = False
        ct.recon.displayImagePictureTitles = True
        retval = recon.recon(ct)

        assert recon_plt_figure.call_count == ct.recon.sliceCount
        assert recon_plt_axis.call_count == ct.recon.sliceCount
        assert recon_plt_imshow.call_count == ct.recon.sliceCount
        assert recon_plt_title.call_count == ct.recon.sliceCount
        assert recon_plt_draw.call_count == ct.recon.sliceCount
        assert recon_plt_savefig.call_count == 0
        assert recon_plt_pause.call_count == 1
        assert recon_plt_close.call_count == 1

    @patch('gecatsim.reconstruction.pyfiles.recon.plt.figure', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.imshow', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.title', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.close', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.draw', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.pause', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.savefig', create=True)
    @patch('gecatsim.reconstruction.pyfiles.recon.plt.axis', create=True)
    def test_recon_save_image_picture_files(self, recon_plt_axis, recon_plt_savefig, recon_plt_pause, recon_plt_draw,
                                          recon_plt_close, recon_plt_title, recon_plt_imshow, recon_plt_figure):
        ct = self.initialize()
        ct.do_Recon = False
        ct.recon.saveImagePictureFiles = True
        ct.waitForKeypress = False
        ct.recon.displayImagePictureTitles = True
        ct.reconImageTitle = 'sample'
        retval = recon.recon(ct)

        assert recon_plt_figure.call_count == ct.recon.sliceCount
        assert recon_plt_axis.call_count == ct.recon.sliceCount
        assert recon_plt_imshow.call_count == ct.recon.sliceCount
        assert recon_plt_title.call_count == ct.recon.sliceCount
        assert recon_plt_draw.call_count == 0
        assert recon_plt_savefig.call_count == ct.recon.sliceCount
        assert recon_plt_pause.call_count == 0
        print(recon_plt_close.call_count)
        assert recon_plt_close.call_count == ct.recon.sliceCount
    def initialize(self):
        ##--------- Initialize
        ct = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")
        ##--------- Make changes to parameters (optional)
        ct.resultsName = "test"
        ct.protocol.viewsPerRotation = 500
        ct.protocol.viewCount = ct.protocol.viewsPerRotation
        ct.protocol.stopViewId = ct.protocol.viewCount - 1
        # ct.protocol.scanTypes = [1, 0, 0]  # flags for airscan, offset scan, phantom scan
        # ct.load_cfg("Protocol_Sample_axial", "Physics_Sample", "Recon_Sample_2d")  # new cfg overrides existing parameters
        ct.protocol.mA = 800
        ct.scanner.detectorRowsPerMod = 4
        ct.scanner.detectorRowCount = ct.scanner.detectorRowsPerMod
        ct.recon.fov = 300.0
        ct.recon.sliceCount = 4  # number of slices to reconstruct
        ct.recon.sliceThickness = 0.568  # reconstruction inter-slice interval (in mm)
        ##--------- Run simulation
        # ct.run_all()  # run the scans defined by protocol.scanTypes
        ##--------- Reconstruction
        ct.do_Recon = True
        return ct


