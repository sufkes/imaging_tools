#!/usr/bin/env python
import sys, os
import argparse

def combineBvecs(path1, path2, outpath):
    with open(path1, 'rb') as fh:
        file1 = fh.readlines()
    with open(path2, 'rb') as fh:
        file2 = fh.readlines()

    lines_combined = []
    for line_num in range(len(file1)):
        lines_combined.append(file1[line_num].rstrip() + ' ' + file2[line_num])

    with open(outpath, 'wb') as fh:
        fh.writelines(lines_combined)
    return

if (__name__ == "__main__"):
    # Create argument parser.
    description = """Input two bval files or two bvec files. 
Concatenates the files and writes them to a new file."""
    parser = argparse.ArgumentParser(description=description)
    
    # Define positional arguments.
    parser.add_argument("path1", type=str, help='path to first bval or bvec file')
    parser.add_argument("path2", type=str, help='path to first bval or bvec file')        
    parser.add_argument("outpath", type=str, help='path to concatenated bval or bvec file')        

    # Print help message if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()
    
    # Parse arguments.
    args = parser.parse_args()

    combineBvecs(args.path1, args.path2, args.outpath)
