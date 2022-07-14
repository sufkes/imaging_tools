#!/usr/bin/env python3

# ISSUES: Currently I am getting errors about misalignment or dimension incongruence.

import os, sys, argparse
import subprocess

from combineBvecs import combineBvecs # for concatenating bval and bvec files

def run_cmd(sys_cmd, verbose=False):
    # One line call to output system command and control debug state.
    if (verbose):
        print(sys_cmd)
    p = subprocess.Popen(sys_cmd, stdout = subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = p.communicate()
    return output, errors

def findFiles(in_dir):
    path_nii = None
    path_bval = None
    path_bvec = None
    for root, dirs, files in os.walk(in_dir):
        for name in files:
            if (name[-4:] == '.nii') or (name[-7:] == '.nii.gz'):
                if (path_nii != None):
                    raise ValueError('Found multiple NIFTI files.')
                path_nii = os.path.join(root, name)
            elif (name[-5:] == '.bval'):
                if (path_bval != None):
                    raise ValueError('Found multiple bval files.')
                path_bval = os.path.join(root, name)
            elif (name[-5:] == '.bvec'):
                if (path_bvec != None):
                    raise ValueError('Found multiple bvec files.')
                path_bvec = os.path.join(root, name)
    return path_nii, path_bval, path_bvec

def combineDTISeries(in_dirs, out_dir):
    # File type and number checks
    if (len(in_dirs) < 2):
        raise Exception("Must specify at least 2 DTI series to combine.")
        return

    # Create out directory
    if (out_dir == None):
        out_dir = os.path.join(os.getcwd(), 'combined')
    if (not os.path.exists(out_dir)):
        print('Creating directory:', out_dir)
        os.mkdir(out_dir)
    elif (not os.path.isdir(out_dir)):
        raise ValueError('out_dir exists, but is not a directory.')
    
    for in_dir in in_dirs:
        # Find paths to files in first directory.
        path_nii, path_bval, path_bvec = findFiles(in_dir)
        if (in_dir == in_dirs[0]): # if first input directory
            ## Copy NIFTI, bval and bvec files to the combined output directory.
            # Set names for the combined files.
            prefix = 'combined' # name to use for combined files
#            if (path_nii[-4] == '.nii'):
#                path_nii_comb = os.path.join(out_dir, prefix+'.nii')
#            else:
#                path_nii_comb = os.path.join(out_dir, prefix+'.nii.gz')
            path_nii_comb = os.path.join(out_dir, prefix+'.nii.gz') # assume fslmerge will use .nii.gz extension
            path_bval_comb = os.path.join(out_dir, prefix+'.bval')
            path_bvec_comb = os.path.join(out_dir, prefix+'.bvec')

            # Copy the files from the first directory.
            cmd_cp_nii = 'cp -n "'+path_nii+'" "'+path_nii_comb+'"'
            out, err = run_cmd(cmd_cp_nii)
            cmd_cp_bval = 'cp -n "'+path_bval+'" "'+path_bval_comb+'"'
            out, err = run_cmd(cmd_cp_bval)
            cmd_cp_bvec = 'cp -n "'+path_bvec+'" "'+path_bvec_comb+'"'
            out, err = run_cmd(cmd_cp_bvec)
        else: # if not first input directory
            ## Combine files.
            # Combine NIFTI files.
            cmd_fslmerge = 'fslmerge -t "'+path_nii_comb+'" "'+path_nii+'" "'+path_nii_comb+'"'
            out, err = run_cmd(cmd_fslmerge)

            # Combine bval and bvec files
            combineBvecs(path_bval_comb, path_bval, path_bval_comb)
            combineBvecs(path_bvec_comb, path_bvec, path_bvec_comb)

if (__name__ == "__main__"):
    # Create argument parser.
    description = """Input two directories. Look for copies of the first directory's 
    images in the second directory. Print a list of all files for which a copy was not 
    found in the second directory. Images for which pixel data matches but header files 
    differ are considered copies"""
    parser = argparse.ArgumentParser(description=description)
    
    # Define positional arguments.
    parser.add_argument("in_dirs", help="paths to directories containing DTI series you wish to combine. Each directory should contain exactly one file ending in '.nii.gz' or '.nii', one file ending in '.bval' and one file ending in '.bvec'.", type=str, nargs='*')
        
    # Define optional arguments.
    parser.add_argument("-o", "--out_dir", action="store", type=str, help="path to directory of combined image and bval/bvec data.")
    
    # Print help message if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()
    
    # Parse arguments.
    args = parser.parse_args()

    combineDTISeries(args.in_dirs, args.out_dir)
