#############################################################################################################
#   Filename:       Caio_Circuit_Solver.py
#   Summary:        Creates an output csv file based on an input .net file describing a cascade circuit
#   Description:    The input file should specify:
#                       - The components of the cirucit
#                       - The frequencies that need to be considered
#                       - The variables the output file should contain
#   Version:        v2.01
#   Date:           28/04/2024
#   Authors:        Caio Pasta (cp2150)
#############################################################################################################

import re
import numpy as np
import pandas as pd
import sys


def split_input_file(input_filename):
    # Opens and reads the input file
    # The function removes comments
    # The lines within each block of the input file are returned

    with open(input_filename, "r") as f:
        contents = f.readlines()

    for i, line in enumerate(contents):
        if "#" in line:
            contents[i] = line.split("#", 1)[0]
        else:
            # Removes the new line character to maintain consistent lines with comments will have them removed
            contents[i] = line.rstrip("\n") 

    circuit_block_start = contents.index("<CIRCUIT>")
    circuit_block_end = contents.index("</CIRCUIT>")

    terms_block_start = contents.index("<TERMS>")
    terms_block_end = contents.index("</TERMS>")

    output_block_start = contents.index("<OUTPUT>")
    output_block_end = contents.index("</OUTPUT>")

    # Incremented to avoid the section indicators
    circuit_block = contents[circuit_block_start+1:circuit_block_end]
    terms_block = contents[terms_block_start+1:terms_block_end]
    output_block = contents[output_block_start+1:output_block_end]

    return circuit_block, terms_block, output_block

def extract_terms(terms_block):
    # Converts Norton Current sources to Thevenin Voltage sources
    # Converts conductance into resistance
    # Calculates the selected frequencies
    # Returns these values 

    # source_data stores the values for conversion to a Thevenin source
    source_data = {"VT":None, "IN":None, "RS":None, "GS":None}
    terms_data = {"VT":None, "RS":None, "RL":None}
    Fstart = None
    Fend = None
    Nfreqs = None
    LFstart = None
    LFend = None
    frequencies = []

    for line in terms_block:
        # Regular expressions are used to efficiently extract the values for each term no matter their format
        match = re.search(r"VT=(-?\d+(\.\d+)?(e[+-]?\d+)?)", line)
        if match:
            source_data["VT"] = float(match.group(1))

        match = re.search(r"IN=(-?\d+(\.\d+)?(e[+-]?\d+)?)", line)
        if match:
            source_data["IN"] = float(match.group(1))

        match = re.search(r"RS=(-?\d+(\.\d+)?(e[+-]?\d+)?)", line)
        if match:
            source_data["RS"] = float(match.group(1))

        match = re.search(r"GS=(-?\d+(\.\d+)?(e[+-]?\d+)?)", line)
        if match:
            source_data["GS"] = float(match.group(1))

        match = re.search(r"RL=(-?\d+(\.\d+)?(e[+-]?\d+)?)", line)
        if match:
            terms_data["RL"] = float(match.group(1))

        match = re.search(r"Fstart=(-?\d+(\.\d+)?(e[+-]?\d+)?)", line)
        if match:
            Fstart = float(match.group(1))

        match = re.search(r"Fend=(-?\d+(\.\d+)?(e[+-]?\d+)?)", line)
        if match:
            Fend = float(match.group(1))

        match = re.search(r"Nfreqs=(-?\d+)", line)
        if match:
            Nfreqs = int(match.group(1))

        match = re.search(r"LFstart=(-?\d+(\.\d+)?(e[+-]?\d+)?)", line)
        if match:
            LFstart = float(match.group(1))

        match = re.search(r"LFend=(-?\d+(\.\d+)?(e[+-]?\d+)?)", line)
        if match:
            LFend = float(match.group(1))

    # Converts the Norton current source into a Thevenin's source
    if source_data["VT"] == None:
        if source_data["RS"] == None:
            source_data["RS"] = 1/source_data["GS"]
            terms_data["VT"] = source_data["IN"] * source_data["RS"]
            terms_data["RS"] = source_data["GS"]
        else:
            terms_data["VT"] = source_data["IN"] * source_data["RS"]
            terms_data["RS"] = 1/source_data["RS"]
    else:
        terms_data["VT"] = source_data["VT"]
        if source_data["RS"] == None:
            terms_data["RS"] = 1/source_data["GS"]
        else:
            terms_data["RS"] = source_data["RS"]

    # Calculates the selected frequencies
    # Whether Fstart or LFstart are found in the input file decides whether the sequence of frequencies is normal or logarithmic
    
    if LFstart:
        frequencies = np.logspace(LFstart, LFend, Nfreqs)
    elif Fstart:
        frequencies = np.linspace(Fstart, Fend, Nfreqs)
    else:
        raise Exception("No frequencies selected")

    return terms_data, frequencies

def create_output_dictionary(output_block):
    # Creates the structure for the dictonary that stores the data for the output file

    output_dictionary = {}

    # 10 is used as the first column is shorter than the rest
    output_dictionary["Freq".rjust(10)] = ["Hz".rjust(10)] 

    for line in output_block:
        properties = re.split(r"\s+", line.strip())
        
        if len(properties) > 1:
            if properties[1][:2] == "dB":
                output_dictionary[("|"+ properties[0] +"|").rjust(11)] = [properties[1].rjust(11)]
                output_dictionary[("/_" + properties[0]).rjust(11)] = ["Rads".rjust(11)]
            else:
                output_dictionary[("Re("+ properties[0] +")").rjust(11)] = [properties[1].rjust(11)]
                output_dictionary[("Im("+ properties[0] +")").rjust(11)] = [properties[1].rjust(11)]
        else:
            # If the property doesn't have a unit specified "L" is used instead
            output_dictionary[("Re("+ properties[0] +")").rjust(11)] = ["L".rjust(11)]
            output_dictionary[("Im("+ properties[0] +")").rjust(11)] = ["L".rjust(11)]

    return output_dictionary

def get_prefix_multiplier(prefix_string):
    # Returns a multiplier based on the prefix provided

    multiplier = 1

    match prefix_string:
        case "p":
            multiplier = 1e-12
        case "n":
            multiplier = 1e-9
        case "u":
            multiplier = 1e-6
        case "m":
            multiplier = 1e-3
        case "k":
            multiplier = 1e3
        case "M":
            multiplier = 1e6
        case "G":
            multiplier = 1e9

    return multiplier

def calculate_ABCD_Matrix(circuit_block, freq):
    # Decodes the circuit block
    # Sorts the components in order from the source to the load
    # Returns the ABCD Matrix of the circuit

    # Converts each line in the Circuit Block to a list
    component_array = []
    for line in circuit_block:
        properties = re.split(r"\s+", line.strip())
        component = [None, None, None, None]
        for p in properties:
            p_segments = p.split("=")

            match p_segments[0]:
                case "n1":
                    component[0] = int(p_segments[1])
                case "n2":
                    component[1] = int(p_segments[1])
                case "R":
                    if len(properties) > 3:
                        prefix_multiplier = get_prefix_multiplier(properties[3].strip())
                    else:
                        prefix_multiplier = 1
                    component[2] = float(p_segments[1]) * prefix_multiplier
                    component[3] = 1/(float(p_segments[1]) * prefix_multiplier)
                case "G":
                    if len(properties) > 3:
                        prefix_multiplier = get_prefix_multiplier(properties[3].strip())
                    else:
                        prefix_multiplier = 1
                    component[2] = 1/(float(p_segments[1]) * prefix_multiplier)
                    component[3] = float(p_segments[1]) * prefix_multiplier
                case "C":
                    if len(properties) > 3:
                        prefix_multiplier = get_prefix_multiplier(properties[3].strip())
                    else:
                        prefix_multiplier = 1
                    # Reversed order to make it easier to read and slightly more efficient
                    component[3] = 2j*np.pi*freq * float(p_segments[1]) * prefix_multiplier
                    component[2] = 1/component[3]
                case "L":
                    if len(properties) > 3:
                        prefix_multiplier = get_prefix_multiplier(properties[3].strip())
                    else:
                        prefix_multiplier = 1
                    component[2] = 2j*np.pi*freq * float(p_segments[1]) * prefix_multiplier
                    component[3] = 1/component[2]

        # If data about a component is missing that component should not be included
        if None not in component:
            component_array.append(component)

    # Sorts the 2D array of components
    component_array.sort(key=lambda x: (x[0], x[1]))

    # Creates an ABCD matrix for the first component
    abcd_Matrix = np.identity(2)
        
    for component in component_array:
        if component[1] == 0:
            temp_matrix = np.array([[1, 0], 
                                    [component[3], 1]])
        else:
            temp_matrix = np.array([[1, component[2]], 
                                    [0, 1]])
        abcd_Matrix = np.matmul(abcd_Matrix, temp_matrix)

    return abcd_Matrix

def calculate_output_values(term_values, abcd_matrix, freq):
    # Returns a dictionary with all the calculated values

    Zin = (abcd_matrix[0][0]*term_values["RL"] + abcd_matrix[0][1])/(abcd_matrix[1][0]*term_values["RL"] + abcd_matrix[1][1])

    Zout = (abcd_matrix[1][1]*term_values["RS"] + abcd_matrix[0][1])/(abcd_matrix[1][0]*term_values["RS"] + abcd_matrix[0][0])
    
    Vin = (Zin/(term_values["RS"] + Zin)) * term_values["VT"]

    Iin = Vin/Zin

    Av = term_values["RL"]/(abcd_matrix[0][0]*term_values["RL"] + abcd_matrix[0][1])

    Ai = 1/(abcd_matrix[1][0]*term_values["RL"] + abcd_matrix[1][1])

    Ap = np.conjugate(Ai) * Av

    Vout = Vin * Av

    Iout = Iin * Ai

    Pin = np.conjugate(Iin) * Vin

    Pout = Pin * Ap

    calculated_values = {"Freq": freq, "Zin": Zin, "Zout": Zout, "Vin": Vin, "Iin": Iin, "Av": Av, "Ai": Ai, "Vout": Vout, "Iout": Iout, "Pin": Pin, "Pout": Pout, "Ap": Ap}

    return calculated_values

def append_to_output_dictionary (calculated_values, output_dictionary):
    # Adds the calculated values to the output dictionary
    # Converting them to the correct format in the process (correct number of significant figures, decibles etc.)
    # Returns the output_dictionary 

    for key in output_dictionary:
        trimmed_key = key.strip()

        if trimmed_key[:2] == "Fr":
            prefix_multiplier = get_prefix_multiplier(output_dictionary[key][0][0])
            corrected_value = calculated_values["Freq"] / prefix_multiplier
            output_dictionary[key].append(("%1.3e" %corrected_value).rjust(10)) # Formats the number to the correct form

        elif trimmed_key[:2] == "Re":
            prefix_multiplier = get_prefix_multiplier(output_dictionary[key][0][0])
            # trimmed_key[3:-1] removes the Re() characters around the variable
            corrected_value = calculated_values[trimmed_key[3:-1]] / prefix_multiplier
            output_dictionary[key].append(("%1.3e" %corrected_value.real).rjust(11))

        elif trimmed_key[:2] == "Im":
            prefix_multiplier = get_prefix_multiplier(output_dictionary[key][0][0])
            corrected_value = calculated_values[trimmed_key[3:-1]] / prefix_multiplier
            output_dictionary[key].append(("%1.3e" %corrected_value.imag).rjust(11))

        elif trimmed_key[:1] == "|":
            prefix_multiplier = get_prefix_multiplier(output_dictionary[key][0][2])
            # If the variable being converted into decibels starts with "V" or "I" or is "Av" or "Ai" 
            # then the equation 20log10(variable_value) needs to be used instead of 10log10(variable_value)
            if (trimmed_key[1] == "V" 
                or trimmed_key[1] == "I" 
                or trimmed_key[1:-1] == "Av" 
                or trimmed_key[1:-1] == "Ai") :
                corrected_value = (20 * np.log10(np.absolute(calculated_values[trimmed_key[1:-1]])/ prefix_multiplier)) 
            else:
                corrected_value = (10 * np.log10(np.absolute(calculated_values[trimmed_key[1:-1]])/ prefix_multiplier)) 
            output_dictionary[key].append(("%1.3e" %corrected_value).rjust(11))

        elif trimmed_key[:2] == "/_":
            prefix_multiplier = get_prefix_multiplier(output_dictionary[key][0][2])
            corrected_value = calculated_values[trimmed_key[2:]] / prefix_multiplier
            output_dictionary[key].append(("%1.3e" %np.angle(corrected_value)).rjust(11))

    return output_dictionary

def main():
    # Abandons the attempt if an error is raised
    try:
        # The input filename is specified by a command line argument
        circuit_block, terms_block, output_block = split_input_file(sys.argv[1])
        terms_data, frequencies = extract_terms(terms_block)
        output_dictionary = create_output_dictionary(output_block)

        # The header dataframe is created with the variables 
        # and then converted to csv
        header_df = pd.DataFrame(output_dictionary)
        header_text = header_df.to_csv(index=False, lineterminator="\n")

        for freq in frequencies:
            abcd_matrix = calculate_ABCD_Matrix(circuit_block, freq)
            calculated_values = calculate_output_values(terms_data, abcd_matrix, freq)
            output_dictionary = append_to_output_dictionary(calculated_values, output_dictionary)

        df = pd.DataFrame(output_dictionary)
        data_text = df.to_csv(index=False, lineterminator=",\n")

        # The previously made header csv is needed to remove the trailing comma
        output_text = header_text.split("\n")[:-1] + data_text.split("\n")[2:-1]

        # The output filename is specified by a command line argument
        with open(sys.argv[2], "wt") as f:
            for line in output_text:
                f.write(line + "\n")

    # If an error is raised the output file should still be created but left empty
    except:
        f = open(sys.argv[2], "wt")
        f.close()

if __name__ == '__main__':
    main()