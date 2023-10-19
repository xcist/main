#import pytest
import gecatsim.pyfiles.CommonTools as c
import numpy as np
import tempfile
import os
from unittest.mock import MagicMock
from unittest.mock import patch,  Mock, call
import numpy.matlib as nm
from io import BytesIO, StringIO

def test_make_col():
    row = 3
    col = 5
    num = 10
    x = np.full((row, col), num)
    x_as_col = c.make_col(x)
    assert len(x_as_col) == row * col


def test_load_C_lib():
    clib = c.load_C_lib()
    assert clib is not None
    print(clib)


def test_path_helper():
    assert c.my_path is not None
    with tempfile.TemporaryDirectory() as tmpdirname:
        with open(os.path.join(tmpdirname, 'temp_phantom.cfg'), 'w') as fp:
            fp.write('Hello world!')
        c.my_path.add_search_path(tmpdirname)
        found_path = c.my_path.find('phantom', 'temp_phantom', '.cfg')
        assert found_path == os.path.join(tmpdirname, 'temp_phantom.cfg')


def test_vector_norm():
    bad_vector = np.zeros([4, 1], dtype=np.single)
    result = c.vectornorm(bad_vector)
    assert result is None

    good_vector = np.ones([3, 1], dtype=np.single)
    result = c.vectornorm(good_vector)
    assert result is not None
    print(result)
    assert abs(result[0][0]-1.7320508) < 1e-5

def test_feval():
    cfg = c.CFG()
    assert len(vars(cfg.det)) == 0
    c.feval('Detector_ThirdgenCurved', cfg);
    assert len(vars(cfg.det)) == 16

    assert cfg.det.nCells is not None
    assert cfg.det.cellCoords  is not None
    assert cfg.det.nSamples  is not None
    assert cfg.det.sampleCoords  is not None
    assert cfg.det.weights is not None
    assert cfg.det.activeArea  is not None
    assert cfg.det.nMod  is not None
    assert cfg.det.modCoords  is not None
    assert cfg.det.uvecs  is not None
    assert cfg.det.vvecs  is not None
    assert cfg.det.totalNumCells  is not None
    assert cfg.det.startIndices  is not None
    assert cfg.det.nModDefs  is not None
    assert cfg.det.modTypes  is not None
    assert cfg.det.width  is not None
    assert cfg.det.height  is not None

def test_feval_invalid_funcname():
    cfg = c.CFG()
    try:
        c.feval('invalid_funcname', cfg);
        assert False
    except:
        assert True

def test_feval_invalid_config():
    cfg = c.emptyCFG()
    try:
        c.feval('invalid_funcname', cfg);
        assert False
    except:
        assert True

def test_path_helper_load():
    #.gecatsim may not be always availablem so not testing the file loading
    pathhelper = c.PathHelper();
    assert len(pathhelper.paths) == 10

    assert pathhelper.paths["main"] is not None
    assert pathhelper.paths["top"] is not None
    assert pathhelper.paths["cfg"].endswith('cfg')
    assert pathhelper.paths["lib"].endswith('lib')
    assert pathhelper.paths["bowtie"].endswith('bowtie')
    assert pathhelper.paths["material"].endswith('material')
    assert pathhelper.paths["phantom"].endswith('phantom')
    assert pathhelper.paths["scatter"].endswith('scatter')
    assert pathhelper.paths["spectrum"].endswith('spectrum')
    assert pathhelper.paths["dose_data"].endswith('dose_data')

def test_find_unavailable_file():
    try:
        found_path = c.my_path.find('phantom', 'temp_phantom', '.cfg')
        assert False
    except:
        assert True

def test_add_dir_to_path():
    tmpdirname = tempfile.TemporaryDirectory().name
    c.PathHelper().add_dir_to_path(tmpdirname)
    assert os.environ["PATH"].find(tmpdirname) != -1
    assert os.environ["PATH"].find(tempfile.TemporaryDirectory().name) == -1

def test_linux_style_path():
    linux_style_path = c.PathHelper().linux_style_path("C:\\Users\\username\\Desktop\\cert.crt");
    assert linux_style_path == "C:/Users/username/Desktop/cert.crt"

def test_cfg():
    cfg = c.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic", "../examples/cfg/Protocol_Sample_axial")
    assert len(vars(cfg)) == 14

    assert cfg.clib is not None
    assert cfg.det is not None
    assert cfg.detNew is not None
    assert cfg.dose is not None
    assert cfg.phantom is not None
    assert cfg.physics is not None
    assert cfg.protocol is not None
    assert cfg.recon is not None
    assert cfg.resultsName is not None
    assert cfg.scanner is not None
    assert cfg.sim is not None
    assert cfg.spec is not None
    assert cfg.src is not None
    assert cfg.srcNew is not None

def test_overlap_no_intersection_x0_lesser():
    x0 = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    x1 = np.array([ 5.0 ,5.5, 6.0, 6.5])
    y0 = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    y1 = np.array([0.0, 0.0, 0.0, 0.0])

    result = c.overlap(x0, y0, x1)

    assert np.array_equal(result, y1)

def test_overlap_no_intersection_x0_higher():
    x0 = np.array([ 5.0 ,5.5, 6.0, 6.5, 7.0, 7.5])
    x1 = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    y0 = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    y1 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    result = c.overlap(x0, y0, x1)

    assert np.array_equal(result, y1)

def test_overlap_no_intersection():
    x0 = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    x1 = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
    y0 = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    y1 = np.array([1.5, 2.0, 2.5, 3.0, 0.0, 0.0])

    result = c.overlap(x0, y0, x1)

    assert np.array_equal(result, y1)

def test_get_vector_boundaries():
    x = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    expected = np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25])

    result = c.get_vector_boundaries(x)

    assert np.array_equal(result, expected)

def test_get_vector_boundaries_single_val():
    x = np.array([1.0])
    expected = np.array([[0.999999], [1.000001]])

    result = c.get_vector_boundaries(x)

    assert np.array_equal(result, expected)

@patch("builtins.open", create=True)
def test_rawread(file_Reader_mock):
    file_mock = BytesIO(b'\x00\x01\x00\x01\x00\x01\x00\x01')
    file_Reader_mock.return_value = file_mock

    c.rawread("samplefile.img", [], 'float')
    assert file_Reader_mock.call_count == 1


@patch("builtins.open", create=True)
def test_rawwrite(file_mock):
    sample_file = CustomStringIO();
    file_mock.return_value = sample_file
    c.rawwrite("sample.txt", "sample_data");

    assert sample_file.getvalue() == "sample_data"

class CustomStringIO(StringIO):
    def close(self):
        pass

