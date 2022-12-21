#!/usr/bin/env python3

import os
import sys
import argparse

import pandas as pd

def main(in_path, ivs, covs, out_dir):
    df = pd.read_csv(in_path)

    # Select only the columns which will be used.
    if covs is None:
        covs = []
    df = df[['filename'] + ivs + covs]

    # Remove rows with blanks.
    df.dropna(axis=0, inplace=True)

    # Save list of files in remaining data.
    file_list_df = df['filename']
    file_list_name = 'file_list.csv'
    file_list_path = os.path.join(out_dir, file_list_name)
    file_list_df.to_csv(file_list_path, index=False, header=False)

    # Save design matrix.
    design_df = df[ivs + covs]
    design_name = 'design.csv'
    design_path = os.path.join(out_dir, design_name)
    design_df.to_csv(design_path, index=False, header=False)

    # Save contrast matrix.
    contrast_df = pd.DataFrame(columns=ivs + covs)
    for ii in range(len(ivs)):
        contrast_df.loc[ii, :] = 0
        contrast_df.iloc[ii, ii] = 1
    contrast_name = 'contrast.csv'
    contrast_path = os.path.join(out_dir, contrast_name)
    contrast_df.to_csv(contrast_path, index=False, header=False)
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Grab specified columns from spreadsheet, remove other columns, remove missing values, save as design matrix for input to FSL PALM; save corresponding contrast file."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
#    parser.add_argument("", help="")
    
    # Define optional arguments.
    parser.add_argument("-i", "--in_path", type=str, help='path to CSV file with data', default='data_sheets/data-master.csv')
    parser.add_argument("-v", '--ivs', help="variables to include in model; contrasts will be included for these variables", nargs='+', type=str)
    parser.add_argument("-c", '--covs', help="covariates to include in model; contrasts will not be included for these variables", nargs='*', type=str)
    parser.add_argument("-o", "--out_dir", type=str, help='output directory')
    
    # Print help if no arguments input.
    #if (len(sys.argv) == 1):
    #    parser.print_help()
    #    sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))
