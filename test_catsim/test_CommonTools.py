import pytest
import catsim.pyfiles.CommonTools as c
import numpy as np
import tempfile
import os
import numpy.matlib as nm


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
