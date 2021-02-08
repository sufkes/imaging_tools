#!/usr/bin/env python

import os
import sys
import argparse

if (__name__ == '__main__'):
    # Create argument parser.
    description = """"""
    parser = argparse.ArgumentParser(description=description)
    
    # Define positional arguments.
#    parser.add_argument("", help="")
    
    # Define optional arguments.
#    parser.add_argument("-", "--", help="")

    # Print help if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args(**var(args))

    # Run main function.
