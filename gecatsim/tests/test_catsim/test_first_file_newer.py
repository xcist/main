import unittest
from unittest.mock import patch
from gecatsim.pyfiles.first_file_newer import first_file_newer, first_date_is_more_recent, monotonic_with_time, numeric_month

class TestFirstFileNewer(unittest.TestCase):

    @patch('gecatsim.pyfiles.first_file_newer.os.path.exists')
    @patch('gecatsim.pyfiles.first_file_newer.os.path.getmtime')
    @patch('gecatsim.pyfiles.first_file_newer.monotonic_with_time')
    def test_first_file_newer_does_not_exist(self, mock_monotonic, mock_getmtime, mock_exists):

        mock_monotonic.side_effect = lambda x: x

        # Test case where first_file does not exist
        mock_exists.side_effect = lambda x: x == 'file2'
        result = first_file_newer('file1', 'file2')
        assert result == False

    @patch('gecatsim.pyfiles.first_file_newer.os.path.exists')
    @patch('gecatsim.pyfiles.first_file_newer.os.path.getmtime')
    @patch('gecatsim.pyfiles.first_file_newer.monotonic_with_time')
    def test_first_file_exists_second_file_does_not(self, mock_monotonic, mock_getmtime, mock_exists):

        mock_monotonic.side_effect = lambda x: x

        # Test case where first_file exists and second_file does not
        mock_exists.side_effect = lambda x: x == 'file1'
        result = first_file_newer('file1', 'file2')
        assert result == True

    @patch('gecatsim.pyfiles.first_file_newer.os.path.exists')
    @patch('gecatsim.pyfiles.first_file_newer.os.path.getmtime')
    @patch('gecatsim.pyfiles.first_file_newer.monotonic_with_time')
    def test_both_files_exist_first_file_newer(self, mock_monotonic, mock_getmtime, mock_exists):

        mock_monotonic.side_effect = lambda x: x

        # Test case where both files exist and first_file is newer
        mock_exists.side_effect = lambda x: True
        mock_getmtime.side_effect = lambda x: 200 if x == 'file1' else 100
        result = first_file_newer('file1', 'file2')
        assert result == True

    @patch('gecatsim.pyfiles.first_file_newer.os.path.exists')
    @patch('gecatsim.pyfiles.first_file_newer.os.path.getmtime')
    @patch('gecatsim.pyfiles.first_file_newer.monotonic_with_time')
    def test_both_files_exist_second_file_newer(self, mock_monotonic, mock_getmtime, mock_exists):

        mock_monotonic.side_effect = lambda x: x

        # Test case where both files exist and second_file is newer
        mock_exists.side_effect = lambda x: True
        mock_getmtime.side_effect = lambda x: 100 if x == 'file1' else 200
        result = first_file_newer('file1', 'file2')
        assert result == False

class TestFirstDateIsMoreRecent(unittest.TestCase):

    @patch('gecatsim.pyfiles.first_file_newer.monotonic_with_time')
    def test_first_date_is_more_recent(self, mock_monotonic):

        mock_monotonic.side_effect = lambda x: x

        # Test case where date1 is more recent than date2
        result = first_date_is_more_recent(200, 100)
        assert result == True

        # Test case where date2 is more recent than date1
        result = first_date_is_more_recent(100, 200)
        assert result == False

class TestMonotonicWithTime(unittest.TestCase):

    def test_monotonic_with_time_valid_date(self):
        # Test case with valid date format
        result = monotonic_with_time('01-Jan-2024 00:00:00')
        assert result == 65052979200

    def test_monotonic_with_time_invalid_date(self):
        # Test case with invalid date format
        with self.assertRaises(ValueError):
            monotonic_with_time('invalid-date')

class TestNumericMonth(unittest.TestCase):

    def test_numeric_month(self):
        # Test case for each month
        assert numeric_month('Jan') == 1
        assert numeric_month('Feb') == 2
        assert numeric_month('Mar') == 3
        assert numeric_month('Apr') == 4
        assert numeric_month('May') == 5
        assert numeric_month('Jun') == 6
        assert numeric_month('Jul') == 7
        assert numeric_month('Aug') == 8
        assert numeric_month('Sep') == 9
        assert numeric_month('Oct') == 10
        assert numeric_month('Nov') == 11
        assert numeric_month('Dec') == 12

if __name__ == '__main__':
    unittest.main()