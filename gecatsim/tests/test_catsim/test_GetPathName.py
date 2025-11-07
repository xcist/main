import unittest
import os
import sys
from gecatsim.pyfiles.GetPathName import GetPathName

class TestGetPathName(unittest.TestCase):

    def test_GetPathName_existing_file(self):
        # Create a temporary file in one of the paths
        temp_dir = sys.path[0]
        temp_file = os.path.join(temp_dir, 'temp_test_file.txt')
        with open(temp_file, 'w') as f:
            f.write('test')

        PathName = GetPathName('temp_test_file.txt')
        self.assertEqual(PathName, temp_file)

        # Clean up
        os.remove(temp_file)

    def test_GetPathName_non_existing_file(self):
        PathName = GetPathName('non_existing_file.txt')
        self.assertEqual(PathName, '')

if __name__ == '__main__':
    unittest.main()