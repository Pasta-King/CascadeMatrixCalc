import unittest
import numpy as np
from Caio_Circuit_Solver import *


# End Testing Class =================================================================
class TestCircuitSolver(unittest.TestCase):

    def test_split_input_file (self):
        test_files = ["test_1", "test_2"]
        expected_outputs = [[["1", "2"], ["5", "6"], ["9", "10", "11"]],
                             [["1", "2"], [], ["7", "8", "9"]]]

        for index, filename in enumerate(test_files):
            circuit_lines, term_lines, output_lines = split_input_file("Personal_Test_Files\\" + filename + ".net")
            self.assertListEqual(circuit_lines, expected_outputs[index][0])
            self.assertListEqual(term_lines, expected_outputs[index][1])
            self.assertListEqual(output_lines, expected_outputs[index][2])

        with self.assertRaises(FileNotFoundError):
            split_input_file("Personal_Test_Files\\non_existant_file.net")

        fail_test_files = ["fail_test_1", "fail_test_2", "fail_test_3"]

        with self.assertRaises(ValueError):
            for filename in fail_test_files:
                split_input_file("Personal_Test_Files\\" + filename + ".net")

    def test_extract_terms(self):
        terms_test_data = ["VT=50", "Fend=50", "RS=1.2e+2", "Fstart=10", "",  "Nfreqs=5", "RL=10"]
        output_dictionary, output_frequencies = extract_terms(terms_test_data)
        expected_dictionary = {'VT': 50.0, 'RS': 120.0, 'RL': 10}
        expected_frequencies = np.array([10., 20., 30., 40., 50.])

        self.assertDictEqual(output_dictionary, expected_dictionary)
        np.testing.assert_array_equal(output_frequencies, expected_frequencies)

        terms_test_data = ["VT=50", "GS=1.2e-2 LFend=3", "LFstart=1",  "Nfreqs=3", "RL=1.2e+02"]
        output_dictionary, output_frequencies = extract_terms(terms_test_data)
        expected_dictionary = {'VT': 50.0, 'RS': 83.3333333, 'RL': 120.}
        expected_frequencies = np.array([10., 100., 1000.])

        for i in output_dictionary:
            self.assertAlmostEqual(output_dictionary[i], expected_dictionary[i])
        np.testing.assert_almost_equal(output_frequencies, expected_frequencies)

        terms_test_data = ["IN=5000e-2", "Nfreqs=5", "RS=50", "Fstart=10",  "Fend=50"]
        output_dictionary, output_frequencies = extract_terms(terms_test_data)
        expected_dictionary = {'VT': 2500.0, 'RS': 0.02, 'RL': None}
        expected_frequencies = np.array([10., 20., 30., 40., 50.])

        for i in output_dictionary:
            self.assertAlmostEqual(output_dictionary[i], expected_dictionary[i])
        np.testing.assert_almost_equal(output_frequencies, expected_frequencies)

        terms_test_data = ["VT=50", "Fend=50", "RS=1.2e+2", "LFstart=10",  "Nfreqs=5"]
        with self.assertRaises(TypeError):
            extract_terms(terms_test_data)
        
        terms_test_data = ["VT=50", "Fend=50", "RS=1.2e+2", "Fstart=10",  "Nfreqs=five"]
        with self.assertRaises(TypeError):
            extract_terms(terms_test_data)

        terms_test_data = ["Fend=50", "RS=1.2e+2", "LFstart=10",  "Nfreqs=5"]
        with self.assertRaises(TypeError):
            extract_terms(terms_test_data)


    def test_create_output_dictionary(self):
        output_test_data = ["Vin", "Vout V", "Iout A", "Pin kW", "Zout mOhms"]
        expected_dict = {"      Freq": ["        Hz"], 
                         "    Re(Vin)": ["          L"], 
                         "    Im(Vin)": ["          L"], 
                         "   Re(Vout)": ["          V"], 
                         "   Im(Vout)": ["          V"], 
                         "   Re(Iout)": ["          A"], 
                         "   Im(Iout)": ["          A"], 
                         "    Re(Pin)": ["         kW"], 
                         "    Im(Pin)": ["         kW"], 
                         "   Re(Zout)": ["      mOhms"],
                         "   Im(Zout)": ["      mOhms"]}

        output_dictionary = create_output_dictionary(output_test_data)

        self.assertDictEqual(output_dictionary, expected_dict)

        output_test_data = ["Vin", "Vout V", "Iout dBA", "Pin kW", "Zout mOhms"]
        expected_dict = {"      Freq": ["        Hz"], 
                         "    Re(Vin)": ["          L"], 
                         "    Im(Vin)": ["          L"], 
                         "   Re(Vout)": ["          V"], 
                         "   Im(Vout)": ["          V"], 
                         "     |Iout|": ["        dBA"], 
                         "     /_Iout": ["       Rads"], 
                         "    Re(Pin)": ["         kW"], 
                         "    Im(Pin)": ["         kW"], 
                         "   Re(Zout)": ["      mOhms"],
                         "   Im(Zout)": ["      mOhms"]}

        output_dictionary = create_output_dictionary(output_test_data)

        self.assertDictEqual(output_dictionary, expected_dict)

    def test_get_prefix_multiplier(self):
        string_inputs = ["m", "V", "dB", "", "kA", "k"]
        expected_outputs = [1e-3, 1, 1, 1, 1, 1e3]

        for idx, x in enumerate(string_inputs):
            output_multiplier = get_prefix_multiplier(x)

            self.assertEqual(output_multiplier, expected_outputs[idx])

    def test_calculate_ABCD_Matrix(self):
        test_circuit_block = ["n1=1 n2=2 R=10",
                              "n1=2 n2=3 C=1.8e-3",
                              "n1=2 n2=0 G=0.5",
                              "n1=3 n2=4 L=3.2 u"]
        expected_output = [[6., 10. -3.51868094j],
                           [0.5, 1. -0.29322341j]]
        
        output_matrix = calculate_ABCD_Matrix(test_circuit_block, 150)

        np.testing.assert_almost_equal(output_matrix, expected_output)

        test_circuit_block = ["n1=1 n2=2 R=10",
                              "n1=2 n2=3 R=0.01 k",
                              "",
                              "n1=2 n2=0 G=0.2 Ohms",
                              "n1=2 n2=0 R=10000 m"]
        expected_output = [[4., 50.],
                           [0.3, 4.]]
        
        output_matrix = calculate_ABCD_Matrix(test_circuit_block, 150)

        np.testing.assert_almost_equal(output_matrix, expected_output)

        with self.assertRaises(ValueError):
            test_circuit_block = ["n1= n2=2 R=10"]
            output_matrix = calculate_ABCD_Matrix(test_circuit_block, 150)
        
        with self.assertRaises(ValueError):
            test_circuit_block = ["n1=1 n2=junk R=10"]
            output_matrix = calculate_ABCD_Matrix(test_circuit_block, 150)

    def test_append_to_output_dictionary(self):
        test_data_dictionary = {"      Freq": ["        Hz"], 
                                "    Re(Vin)": ["          V"], 
                                "    Im(Vin)": ["          V"],
                                "   Re(Zout)": ["      Ohms"],
                                "   Im(Zout)": ["      Ohms"]}
        test_calculated_data = {"Freq":150, "Vin": 25, "Zout": 5e+4 - 200j}

        model_data_dictionary = {"      Freq": ["        Hz", " 1.500e+02"], 
                                "    Re(Vin)": ["          V", "  2.500e+01"], 
                                "    Im(Vin)": ["          V", "  0.000e+00"],
                                "   Re(Zout)": ["      Ohms", "  5.000e+04"], 
                                "   Im(Zout)": ["      Ohms", " -2.000e+02"]}

        output_data_dictionary = append_to_output_dictionary(test_calculated_data, test_data_dictionary)
        self.assertDictEqual(output_data_dictionary, model_data_dictionary)

        test_data_dictionary = {"      Freq": ["        Hz"], 
                                "    Re(Vin)": ["          V"], 
                                "    Im(Vin)": ["          V"],
                                "     |Zout|": ["     dBOhms"],
                                "     /_Zout": ["       Rads"]}
        test_calculated_data = {"Freq":150, "Vin": 25, "Zout": 5e+4 - 200j}

        model_data_dictionary = {"      Freq": ["        Hz", " 1.500e+02"], 
                                "    Re(Vin)": ["          V", "  2.500e+01"], 
                                "    Im(Vin)": ["          V", "  0.000e+00"],
                                "     |Zout|": ["     dBOhms", "  4.699e+01"],
                                "     /_Zout": ["       Rads", " -4.000e-03"]}

        output_data_dictionary = append_to_output_dictionary(test_calculated_data, test_data_dictionary)
        self.assertDictEqual(output_data_dictionary, model_data_dictionary)

# End Testing Class =================================================================

if __name__ == '__main__':
    unittest.main()