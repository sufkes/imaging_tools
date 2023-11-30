#!/usr/bin/env python3

import os
import sys
import argparse

import pandas as pd

def main(request_path, census_path, request_postal_code_column_name, census_postal_code_column_name, fields_to_query, out_path, census_field_column, census_value_column):
    request_df = pd.read_csv(request_path)
    census_df = pd.read_csv(census_path)

    # Remove all of the unnecessary rows from the giant census file.
    census_df = census_df.loc[census_df[census_field_column].isin(fields_to_query), :]
    
    # Loop over the requested postal codes and get the corresponding census fields. 
    for index in request_df.index:
        postal_code = request_df.loc[index, request_postal_code_column_name]
        for field in fields_to_query:
            value = census_df.loc[(census_df[census_postal_code_column_name]==postal_code) & (census_df[census_field_column]==field), census_value_column]
            if len(value) == 0:
                continue
            elif len(value) == 1:
                value = value.values[0]
                request_df.loc[index, field] = value
            else:
                print(f'Multiple values returned for field: {field}')
                request_df.loc[index,field] = 'Error: Multiple values found'

    # Save requested postal codes and corresponding census fields to CSV.
    request_df.to_csv(out_path, index=False)
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Get data from Canadian census corresponding to forward sortation areas (FSA; the first 3 digits of postal code). Input a CSV with a column of FSAs, and a census file broken down by FSA. The census data was downloaded from here: https://www12.statcan.gc.ca/census-recensement/2016/dp-pd/prof/details/download-telecharger/comp/page_dl-tc.cfm?Lang=E -> Forward sortation areas (FSAs) -> CSV"""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
#    parser.add_argument("", help="")
    
    # Define optional arguments.
    parser.add_argument("--request_path", help="path to CSV with a column of 3 digit postal codes for which we want to extract data from the census file", default="postal_codes_split.csv")
    parser.add_argument("--census_path", help="path to census file", default="98-401-X2016046_English_CSV_data.csv")
    parser.add_argument("--request_postal_code_column_name", help="name of column in the requested postal code file containing the first 3 digits of postal codes", default="postal_code")
    parser.add_argument("--census_postal_code_column_name", help="name of column in the census file containing the first 3 digits of postal codes", default="GEO_CODE (POR)")
    parser.add_argument("-f", "--fields_to_query", help="datapoints to return for each requested postal code, as they appear in column 9 of the census file.", nargs="+", default=["Median total income of households in 2015 ($)"], metavar=('FIELD_1', 'FIELD_2'))
    parser.add_argument("-o", "--out_path", help="path to output CSV file", default="postal_codes_split-with_census_data.csv")
    parser.add_argument("--census_field_column", help="name of column in the census file that contains the names of the data fields", default="DIM: Profile of Forward Sortation Areas (2247)")
    parser.add_argument("--census_value_column", help="name of column in the census file that contains the values of the data fields", default="Dim: Sex (3): Member ID: [1]: Total - Sex")
                    
    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))


