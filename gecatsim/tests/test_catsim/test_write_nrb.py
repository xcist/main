import unittest
import numpy as np
from gecatsim.pyfiles.write_nrb import nrbwrite

class TestNrbWrite(unittest.TestCase):

    def test_nrbwrite(cls):
        filename = 'test_output.txt'
        xp = [np.array([[0, 1], [2, 3]])]
        yp = [np.array([[4, 5], [6, 7]])]
        zp = [np.array([[8, 9], [10, 11]])]
        uvecs = [np.array([0.0, 0.5, 1.0])]
        vvecs = [np.array([0.0, 0.5, 1.0])]
        mat_ind = [1]
        obj_names = ['TestObject']

        nrbwrite(filename, xp, yp, zp, uvecs, vvecs, mat_ind, obj_names)

        with open(filename, 'r') as f:
            content = f.readlines()

        expected = [
            '\n',
            'TestObject\n',
            '1\n',
            '2 :M\n',
            '2 :N\n',
            'U Knot Vector\n',
            '0.000000\n',
            '0.500000\n',
            '1.000000\n',
            'V Knot Vector\n',
            '0.000000\n',
            '0.500000\n',
            '1.000000\n',
            'Control Points\n',
            '0.000000 4.000000 8.000000\n',
            '2.000000 6.000000 10.000000\n',
            '1.000000 5.000000 9.000000\n',
            '3.000000 7.000000 11.000000\n'
        ]

        cls.assertEqual(content, expected)

if __name__ == '__main__':
    unittest.main()
