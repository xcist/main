import unittest
from unittest.mock import patch, mock_open
import numpy as np
import os
from gecatsim.pyfiles.generate_subtraction_material import generate_subtraction_material

mt1 = "material1.txt"
mt2 = "material2.txt"
mt_new = "material_new.txt"

@patch("builtins.open", new_callable=mock_open)
@patch("gecatsim.pyfiles.ReadMaterialFile")
@patch("gecatsim.pyfiles.GetMu")
def test_generate_subtraction_material(mock_get_mu, mock_read_material_file, mock_open):
    # Mock the ReadMaterialFile function to return the correct number of values
    mock_read_material_file.side_effect = [
        (2, 1.0, [1, 2], [0.5, 0.5]),
        (2, 1.0, [1, 3], [0.3, 0.7])
    ]

    # Mock the GetMu function
    mock_get_mu.side_effect = [
        np.array([0.1, 0.2, 0.3]),
        np.array([0.05, 0.15, 0.25]),
        np.array([0.05, 0.05, 0.05])
    ]

    # Mock the open function to handle file reading and writing
    m_open = mock_open()
    m_open.side_effect = [
        mock_open(read_data="2\n1.0\n1 0.5\n2 0.5\n").return_value,
        mock_open(read_data="2\n1.0\n1 0.3\n3 0.7\n").return_value,
        mock_open().return_value
    ]

    with patch("builtins.open", m_open):
        AtomicNumbers_new, NormalizedMassFractions_new, Density_new = generate_subtraction_material(mt1, mt2, mt_new)

        assert AtomicNumbers_new is not None
        assert NormalizedMassFractions_new is not None
        assert Density_new is not None

if __name__ == '__main__':
    unittest.main()