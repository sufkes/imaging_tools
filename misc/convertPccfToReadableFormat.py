#!/usr/bin/env python

import os
import sys
import argparse

import numpy as np
import pandas as pd

def main(pccf_path, ref_path, out_path, postal_code_request_path, postal_code_request_header, single_link_indicator_rows_only):
    df_ref = pd.read_csv(ref_path, encoding='latin')
    df_out = pd.DataFrame(columns=df_ref['Attribute_name'])

    # Open the PCCF file as an array.
    with open(pccf_path, 'r', encoding='latin') as handle:
        # Check number of characters in first line (assume same length for all lines).
        line_length = len(handle.readline().strip('\n'))
        array_dtype = f'<U{line_length}'

        # Check that number of characters in line matches the expected length based on the reference table.
        expected_line_length = df_ref['Size'].sum()
        if line_length != expected_line_length:
            raise Exception(f'Lines containing exactly {expected_line_length} characters were expected based on the input PCCF Reference Guide, but first line of PCCF contains {line_length} characters.')
        
        pccf_array = np.fromiter((line for line in handle), dtype=array_dtype)
        
    # Extract variables from each line in the PCCF.
    for row in df_ref.index:
        attribute_name = df_ref.loc[row, 'Attribute_name']
        position = df_ref.loc[row, 'Position']
        size = df_ref.loc[row, 'Size']

        start = position - 1
        stop = start + size
        df_out[attribute_name] = [line[start:stop] for line in pccf_array]

    # Convert datatypes as necessary.
    for column in df_out.columns:
        datatype = df_ref.loc[df_ref['Attribute_name']==column, 'Type'].values[0]
        if datatype == 'C':
            pass # Default is already string.
        elif datatype == 'N':
            df_out[column] = df_out[column].astype(float)

    # Remove rows not marked with the Single Link Indicator (SLI).
    if single_link_indicator_rows_only:
        df_out = df_out[df_out['SLI']=='1']

    # Get the requested postal codes from the complete dataframe.
    if not postal_code_request_path is None:
        if postal_code_request_header:
            header = 0
        else:
            header = None
        df_req = pd.read_csv(postal_code_request_path, header=header, dtype=str)
        df_req.columns = ['Postal code'] # rename the only column in this file.
        df_out = pd.merge(left=df_req, right=df_out, on='Postal code', how='left')
    
    df_out.to_excel(out_path, index=False)
    #df_out.to_csv(out_path, index=False)
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Convert Canadian Postal Code Conversion File (PCCF) to a readable spreadsheet with column headers, optionally excluding rows not in a specified list of postal codes, and optionally returning only 'single link indicator' rows for postal codes with multiple rows in the PCCF."""
    class formatter_class(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(description=description, formatter_class=formatter_class)
    
    # Define positional arguments.
    parser.add_argument("pccf_path", help="path to Postal Code Conversion File (PCCF) text file.")
    parser.add_argument("out_path", help="output path for readable version of PCCF (XLSX)")
    
    
    # Define optional arguments.
    parser.add_argument("--ref_path", help="path to CSV version of Table 4.1 - Record layout and data description from the PCCF Reference Guide. This table specifies the Position, Size (number of characters), Type (C=character, N=Number?), Attribute name, and Description for each variable in the PCCF.", default=os.path.join(os.path.dirname(__file__), 'PCCF_202212-eng-table4.1_PCCF_record_layout.csv'))
    parser.add_argument("--postal_code_request_path", help="path to list of postal codes to be included in output file; should contain one postal code per row.")
    parser.add_argument("--postal_code_request_header", action="store_true", help="use this flag if the postal code request path has a header line which should be ignored")
    parser.add_argument("-s", "--single_link_indicator_rows_only", action="store_true", help="only include records marked as the single link indicator (i.e. rows with SLI=1). This ensures a one-to-one mapping between postal code and geographic region.")

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))
