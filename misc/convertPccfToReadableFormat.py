#!/usr/bin/env python

import os
import sys
import argparse

import numpy as np
import pandas as pd

def main(pccf_path, out_path, ref_path, additional_pccf_paths, single_link_indicator_rows_only, remove_duplicate_postal_codes, postal_code_request_path, postal_code_request_header):
    df_ref = pd.read_csv(ref_path, encoding='latin')
    df_out = pd.DataFrame(columns=df_ref['Attribute_name'])

    # Open the PCCFs file as arrays and concatenate them.
    print('Opening PCCFs')
    if additional_pccf_paths is None:
        pccf_paths = [pccf_path]
    else:
        pccf_paths = [pccf_path] + additional_pccf_paths
    pccf_array = None
    for this_pccf_path in pccf_paths:
        with open(this_pccf_path, 'r', encoding='latin') as handle:
            # Check number of characters in first line (assume same length for all lines).
            line_length = len(handle.readline().strip('\n'))
            array_dtype = f'<U{line_length}'
            
            # Check that number of characters in line matches the expected length based on the reference table.
            expected_line_length = df_ref['Size'].sum()
            if line_length != expected_line_length:
                raise Exception(f'Lines containing exactly {expected_line_length} characters were expected based on the input PCCF Reference Guide, but first line of PCCF contains {line_length} characters.')

            this_pccf_array = np.fromiter((line for line in handle), dtype=array_dtype)
            if pccf_array is None:
                pccf_array = this_pccf_array.copy()
            else:
                pccf_array = np.concatenate((pccf_array, this_pccf_array), axis=0)

    # Extract variables from each line in the PCCF.
    print('Parsing PCCFs')
    for row in df_ref.index:
        attribute_name = df_ref.loc[row, 'Attribute_name']
        position = df_ref.loc[row, 'Position']
        size = df_ref.loc[row, 'Size']

        start = position - 1
        stop = start + size
        df_out[attribute_name] = [line[start:stop] for line in pccf_array]
    postal_code_variable_name = df_out.columns[0] # Get name of postal code variable.
    print(f"Warning: Assuming the first variable in the PCCF, named '{postal_code_variable_name}' represents Postal Code.")

    # Convert datatypes as necessary.
    print('Converting datatypes.')
    for column in df_out.columns:
        datatype = df_ref.loc[df_ref['Attribute_name']==column, 'Type'].values[0]
        if datatype == 'C':
            pass # Default is already string.
        elif datatype == 'N':
            df_out[column] = df_out[column].astype(float)

    # Remove rows not marked with the Single Link Indicator (SLI).
    if single_link_indicator_rows_only:
        print('Removing non-Single-Link-Indicator rows')
        df_out = df_out[df_out['SLI']=='1']

    # Remove rows with duplicate postal codes.
    if remove_duplicate_postal_codes:
        print('Removing duplicate postal codes')
        df_out = df_out[~df_out[postal_code_variable_name].duplicated()]
        
    # Get the requested postal codes from the complete dataframe.
    if not postal_code_request_path is None:
        print('Removing postal codes not in request list.')
        if postal_code_request_header:
            header = 0
        else:
            header = None
        df_req = pd.read_csv(postal_code_request_path, header=header, dtype=str)
        df_req.columns = [postal_code_variable_name] # rename the only column in this file.
        df_out = pd.merge(left=df_req, right=df_out, on=postal_code_variable_name, how='left')

    # Save processed PCCF to a spreadsheet.
    print('Writing to XLSX')
    if len(df_out) >= 1048576: # if too big to save to a single XLSX file, split into pieces.
        max_rows = 1000000
        num_out_files = int(np.ceil(len(df_out)/max_rows))
        print(f"Warning: Output data too large for a single XLSX file; splitting into {num_out_files} files.")
        for i in range(num_out_files):
            start = i*max_rows
            stop = min((i+1)*max_rows, len(df_out))
            this_out_path = out_path[:-5] + f"-{i+1}of{num_out_files}.xlsx"
            this_df_out = df_out.iloc[start:stop, :]
            this_df_out.to_excel(this_out_path, index=False)
    else:    
        df_out.to_excel(out_path, index=False)

    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Convert Canadian Postal Code Conversion File (PCCF) to a readable spreadsheet with column headers, optionally excluding rows not in a specified list of postal codes, and optionally returning only 'single link indicator' rows for postal codes with multiple rows in the PCCF."""
    class formatter_class(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(description=description, formatter_class=formatter_class)
    
    # Define positional arguments.
    parser.add_argument("pccf_path", help="path to Postal Code Conversion File (PCCF) text file.")
    parser.add_argument("out_path", help="output path for readable version of PCCF (must end with .xlsx)")
    
    
    # Define optional arguments.
    parser.add_argument("-r", "--ref_path", help="path to CSV version of Table 4.1 - Record layout and data description from the PCCF Reference Guide. This table specifies the Position, Size (number of characters), Type (C=character, N=Number?), Attribute name, and Description for each variable in the PCCF. By default, use the PCCF Reference Guide released 2022 December 15.", default=os.path.join(os.path.dirname(__file__), 'PCCF_202212-eng-table4.1_PCCF_record_layout.csv'))
    parser.add_argument("-a", "--additional_pccf_paths", help="paths to additional PCCF files (e.g. for retired postal codes) to be appended to the first PCCF. All input PCCFS should have the same format (e.g. mixing PCCFS from 2001 and 2021 will cause an error due to different format).", nargs='+', type=str)
    parser.add_argument("-s", "--single_link_indicator_rows_only", action="store_true", help="only include records marked as the single link indicator (i.e. rows with SLI=1). This ensures a one-to-one mapping between postal code and geographic region.")
    parser.add_argument("-d", "--remove_duplicate_postal_codes", help="remove rows with duplicate postal codes; only the first row in which given postal code appears will be kept. These deletions are made after removing rows with SLI=0 if the --single_link_indicator_rows_only was used.", action='store_true')
    parser.add_argument("-p", "--postal_code_request_path", help="path to list of postal codes to be included in output file; should contain one postal code per row.")
    parser.add_argument("-t", "--postal_code_request_header", action="store_true", help="use this flag if the postal code request path has a header line which should be ignored")

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))
