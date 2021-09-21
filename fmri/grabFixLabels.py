#!/usr/bin/env python3


# E.g. lines in the CSV file made by Julia:
#Subject,MS040027_V02,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#IC,,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,,,,,,,,,,,,,,,,,,,,,,,,,,,
#Score (0-5),,3,3,5,5,5,5,5,5,3,2,0,1,3,0,3,,,,,,,,,,,,,,,,,,,,,,,,,,,

import os
import sys
import argparse

def main(in_path, out_dir, threshold):
    with open(in_path, 'r') as fin:
        lines = fin.readlines()

    for n, line in enumerate(lines):
        if 'Subject' in line:
            subject = line.split(',')[1]

            if subject.strip() == '':
                continue
            
            ics = lines[n+1].split(',')
            scores = lines[n+2].split(',')

            out = '['
            for ii, ic in enumerate(ics):
                try:
                    ic_num = int(ic) # if string can be converted to int without value error, it is a valid IC column.
                    ic_score = int(scores[ii])

                    if ic_score <= threshold:
                        out += str(ic_num)+','
                except ValueError:
                    continue
            out = out.rstrip(',') # remove trailing comma
            out += ']\n'

            print(subject)
            print(out)
            
            # Write output file
            out_name = subject+'.txt'
            out_path = os.path.join(out_dir, out_name)
            with open(out_path, 'w') as fout:
                fout.write(out)
        
in_path = str(sys.argv[1])        
if (__name__ == '__main__'):
    # Create argument parser.
    description = """"""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("in_path", type=str, help="path to CSV file with ratings.")
    
    # Define optional arguments.
    parser.add_argument("-t", "--threshold", help="maximum score to be rated 'noise'", type=int, default=2)
    parser.add_argument("-o", "--out_dir", help="directory to save label text files to", type=str, default=os.getcwd())

    # Print help if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))
