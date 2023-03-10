#!/usr/bin/env python3

import os
import sys
import argparse
from collections import OrderedDict

def main(in_path, out_csv, postal_codes_csv, postal_codes_csv_has_no_header, pc_start, pc_end, da_start, da_end, fill_value_for_unmatched_postal_code):
    ## Read in the conversion file.
    with open(in_path, 'r', encoding='latin') as handle:
        lines = handle.readlines()

    ## Create dictionary mapping from postal code to DAuid
    d = OrderedDict()
    for line in lines:
        postal_code = line[pc_start-1:pc_end]
        dauid = line[da_start-1:da_end]
        d[postal_code] = dauid

    ## Write to CSV file if specified.
    with open(out_csv, 'w') as handle:
        # Add values.
        if postal_codes_csv is None:
            ## Convert all postal codes in the conversion file.
            # Add header.
            header = 'postal_code,dauid\n'
            handle.write(header)
            
            # Convert postal codes to DAuid.
            for postal_code, dauid in d.items():
                handle.write(f'{postal_code},{dauid}\n')
        else:    
            ## Convert the requested postal codes.
            # Open requested postal codes CSV.
            with open(postal_codes_csv, 'r') as request_handle:
                request_lines = request_handle.readlines()
                
            # Add header.
            if not postal_codes_csv_has_no_header:
                header = request_lines[0].strip().split(',')[0] # Assume the requested postal code CSV contains only a single column.
                header += ',dauid\n'
            else:
                header = 'postal_code,dauid\n'
            handle.write(header)
            
            # Convert postal codes to DAuid.
            start_row = 0 if postal_codes_csv_has_no_header else 1
            for line in request_lines[start_row:]:
                postal_code = line.strip().split(',')[0]
                try:
                    dauid = d[postal_code]
                except KeyError:
                    dauid = fill_value_for_unmatched_postal_code
                handle.write(f'{postal_code},{dauid}\n')
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Convert postal code to DAuid"""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("in_path", help="postal code conversion text file")
    parser.add_argument("out_csv", help="path to output CSV file with columns 'postal_code,dauid'", type=str)
    
    # Define optional arguments.
    parser.add_argument("-p", "--postal_codes_csv", type=str, help="path to CSV file containing a single column of postal codes; a column of corresponding DAuids will be added. If this option is not used, the output file will contain every postal code from the conversion text file.")
    parser.add_argument("-n", "--postal_codes_csv_has_no_header", action='store_true', help="use this if the postal_codes_csv file does not have a header, to ensure the first row is also processed.")
    parser.add_argument("--pc_start", type=int, help="start column of postal code (start from 1)", default=1)
    parser.add_argument("--pc_end", type=int, help="start column of postal code (start from 1)", default=6)
    parser.add_argument("--da_start", type=int, help="start column of DAuid (start from 1)", default=126)
    parser.add_argument("--da_end", type=int, help="start column of DAuid (start from 1)", default=133)
    parser.add_argument("--fill_value_for_unmatched_postal_code", type=str, default='', help='DAuid value to use if requested postal code is not found in the conversion file.')
    
    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))
