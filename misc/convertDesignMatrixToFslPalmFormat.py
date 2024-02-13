#!/usr/bin/env python3

import os
import shutil
import sys
import argparse
import warnings

import pandas as pd
import nibabel as nib
import numpy as np

def main(design_matrix_paths, contrast_matrix_paths, contrast_has_row_column_names, out_image_path, file_path_column):
    # Open the design matrix.
    design_df = pd.read_csv(design_matrix_paths[0])
    
    # Generate the 4D NIFTI file.
    if not out_image_path is None:
        if file_path_column is None:
            msg = "If --out_image_path is specified, you must also specify --file_path_column"
            raise Exception(msg)

        file_path_list = design_df[file_path_column].tolist()
        num_images = len(file_path_list)
        first_image = True
        for index, file_path in enumerate(file_path_list):
            subject_nii = nib.load(file_path)
            subject_array = subject_nii.get_fdata()

            # Initialize the 4D image using the first input image's affine and header.
            if first_image:
                combined_affine = subject_nii.affine
                combined_header = subject_nii.header
                combined_shape = (subject_array.shape) + (num_images,) # dimensions 1-3 should match all of the input image dimensions; dimension 4 is the number of input.

                combined_array = np.zeros(shape=combined_shape)
                
                first_image = False

            # Check that dimensions match across inputs.
            if (subject_array.shape != combined_shape[:3]):
                raise Exception(f"Input images have different shapes. First image with a different shape is: {file_path}")
            
            # Check that affines match across inputs.
            if (combined_affine != subject_nii.affine).any():
                warnings.warn(f"Input image may be in a different space than others: {file_path}")
                
            # Add the subject's image to the combined image.
            combined_array[:, :, :, index] = subject_array
        # Save the combined array. 
        combined_nii = nib.nifti1.Nifti1Image(combined_array, affine=combined_affine, header=combined_header)
        nib.save(combined_nii, out_image_path)

    ## Process the design matrix for input to PALM.
    # Remove the subject ID and file path columns.
    design_df.drop(labels=design_df.columns[0], axis=1, inplace=True)
    if file_path_column is not None:
        design_df.drop(labels=file_path_column, axis=1, inplace=True)

    # Save with PALM formatting requirements.
    design_df.to_csv(design_matrix_paths[1], header=False, index=False, encoding='ascii', lineterminator='\n')

    # Save copy of original design matrix for future reference.
    if design_matrix_paths[1].endswith('.csv'):
        copy_of_original_path = f'{design_matrix_paths[1][:-4]}-with_header_before_modification_for_palm.csv'
    else:
        copy_of_original_path = f'{design_matrix_paths[1]}-with_header_before_modification_for_palm.csv'
    shutil.copy2(design_matrix_paths[0], copy_of_original_path)

    ## Process the contrast matrix for input to PALM.
    if contrast_matrix_paths is not None:
        for contrast_matrix_path_pair in contrast_matrix_paths[:]: # there could be more than one contrast matrix to convert (e.g. a t contrast file, and an f contrast file).
            # Open.
            if contrast_has_row_column_names:
                header = 0
            else:
                header = None
            contrast_df = pd.read_csv(contrast_matrix_path_pair[0], header=header)
        
            # Remove row names if necessary.
            if contrast_has_row_column_names:
                contrast_df.drop(labels=contrast_df.columns[0], axis=1, inplace=True)
        
            # Save with PALM formatting requirements.
            contrast_df.to_csv(contrast_matrix_path_pair[1], header=False, index=False, encoding='ascii', lineterminator='\n')
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Convert human-readable design and contrast matrix CSVs into the format required by FSL PALM. Optionally, generate a 4D NIFTI file for input to PALM, using file paths listed in the input design matrix."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
#    parser.add_argument("", help="")
    
    # Define optional arguments.
    parser.add_argument("-d", "--design_matrix_paths", type=str, help="Paths to input design matrix, and the processed design matrix to be used in PALM, both in CSV format. The input should be one subject per row. The subject ID should appear in the first column. Every column other than the first will be treated as variables in the model. Categorical variables should be encoded as 0s and 1s, with a separate column for each group (e.g. sex should be encoded in two columns: one for male and one for female).", nargs=2, metavar=('IN_DESIGN_PATH','OUT_DESIGN_PATH'), required=True)
    parser.add_argument("-c", "--contrast_matrix_paths", type=str, help="Paths to the input contrast matrix, and the processed contrast matrix to be used in PALM, both in CSV format. The input should be one row per test, and one column per variable. This flag can be used multiple times to process multiple contrast files (e.g. a t contrast and an f contrast).", nargs=2, metavar=('IN_CONTRAST_PATH','OUT_CONTRAST_PATH'), action='append')

    #parser.add_argument("--no_names_in_design", action="store_true", help="input design matrix does not have column (variable) names or row (subject) names; do not remove the first row or column") # NOT IMPLEMENTED
    parser.add_argument("--contrast_has_row_column_names", action="store_true", help="Input contrast matrix has column and row names. These will be removed in the processed file.")

    parser.add_argument("-o", "--out_image_path", type=str, help="Generate a 4D NIFTI file with subjects ordered as they appear in the design matrix; the paths to the 3D NIFTI files should appear in design matrix under a specified column. To use this option specify the path of the output 4D NIFTI file here, and specify the name of the column which stores the input file paths in the --file_path_column argument. The column storing the input file paths will be removed from the design matrix. All input images should be in the same space, as required by FSL PALM. For PALM to run, it may be necessary to use extension '.nii' instead of the compressed '.nii.gz' format.", metavar=("OUT_IMAGE_PATH",))
    parser.add_argument("-f", "--file_path_column", type=str, help="If generating a 4D image containing all subjects, specify the name of the column in the design matrix which stores the individual subject image file paths.")

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))
