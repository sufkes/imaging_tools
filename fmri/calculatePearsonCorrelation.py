#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import pandas as pd

def main(time_series_path, out_dir, no_csv=True, scan_id=None, csv_prefix='all_pearson_r', csv_index_col='scan_id'):
    # Load the mean time series.
    mts = np.loadtxt(time_series_path)
    r_raw = np.corrcoef(mts, rowvar=False)

    # numpy.corrcoef may return a non-symmetric result. Force it to be symmetric.
    r_raw = (r_raw + r_raw.T)/2
    if (r_raw != r_raw.T).any():
        print('Warning: Correlation coefficient matrix is not symmetric, despite attempt to force it to be symmetric.')

    r_abs = np.abs(r_raw)
    
    r_pos = r_raw.copy()
    r_pos[r_raw<0] = 0

    # Set the name of the scan.
    if scan_id is None:
        scan_id = os.path.basename(time_series_path)[:12]

    # Create output directory.
    os.makedirs(out_dir, exist_ok=True)
        
    # Save the correlation matrices.
    np.save(os.path.join(out_dir, scan_id+'-raw.npy'), r_raw)
    np.save(os.path.join(out_dir, scan_id+'-abs.npy'), r_abs)
    np.save(os.path.join(out_dir, scan_id+'-pos.npy'), r_pos)

    # Save a CSV file containing the results for all scans.
    if not no_csv:
        csv_path_raw = os.path.join(out_dir, csv_prefix+'-raw.csv')
        csv_path_abs = os.path.join(out_dir, csv_prefix+'-abs.csv')
        csv_path_pos = os.path.join(out_dir, csv_prefix+'-pos.csv')
        if os.path.isfile(csv_path_raw):
            df_raw = pd.read_csv(csv_path_raw, index_col=csv_index_col)
        else:
            df_raw = pd.DataFrame()
            df_raw.index.name = csv_index_col
            
        num_rois = mts.shape[1]

        # Add in the correlation values.
        for roi_1 in range(1, num_rois+1):
            for roi_2 in range(roi_1+1, num_rois+1):
                col_name = f'r_raw_{roi_1}_{roi_2}'
                df_raw.loc[scan_id, col_name] = r_raw[roi_1-1, roi_2-1]

        # Create the absolute value and positive variations of the correlation dataframes.
        df_abs = df_raw.copy()
        df_abs.rename(columns={col:col.replace('raw', 'abs') for col in df_raw.columns}, inplace=True)
        df_abs = df_abs.abs()

        df_pos = df_raw.copy()
        df_pos.rename(columns={col:col.replace('raw', 'pos') for col in df_raw.columns}, inplace=True)
        df_pos[df_pos<0] = 0

        # Save the dataframes as CSV files.
        df_raw.to_csv(csv_path_raw, index=True)
        df_abs.to_csv(csv_path_abs, index=True)
        df_pos.to_csv(csv_path_pos, index=True)
    
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """"""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("time_series_path", help="path to time series text file. The file should be space delimited, time increasing with row number, and label increasing with column number, as output by the FSL command fslmeants with the --label argument specified.")
    parser.add_argument("out_dir", help="path to directory in which to save the correlation matrices and summary CSV files")
    
    # Define optional arguments.
    parser.add_argument("-c", "--no_csv", action='store_true', help="do not add the correlations to a CSV containing data for all scans")
    parser.add_argument("-i", "--scan_id", type=str, help="unique name of scan to be used in the CSV and filenames. By default, use the first 12 characters of the time series file name.")

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))
