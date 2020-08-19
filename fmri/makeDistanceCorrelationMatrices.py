#!/usr/bin/env python

import os, sys
import argparse

if (__name__ == '__main__'):
    # Create argument parser
    description = """Description of function"""
    epilog = '' # """Text to follow argument explantion """
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
#    parser.add_argument("in_path", help="path to DICOM file to be anonymized, or directory containing DICOM files.")
    
    # Define optional arguments.
#    parser.add_argument("-n", "--name", help="Subject ID to set PatientName tags to", type=str)

    # Parse arguments.
    args = parser.parse_args()


