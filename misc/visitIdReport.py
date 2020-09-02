#!/usr/bin/env python

import os, sys
import argparse

import pandas
import glob


def visitIdReport(in_dir, out_path, ignore_symlinks=False, replace_visit_prefix=False):
    report_df = pandas.DataFrame()
    dir_list = glob.glob(in_dir+'/*')
    dir_list.sort()

    for dd in dir_list:
        if (dd.endswith('.txt')):
            continue
        if (not os.path.isdir(dd)):
            if ignore_symlinks and (os.path.islink(dd)):
                print "Ignoring link: ", dd
                continue
            else:
                print "Including non-directory:", dd

        dd = os.path.basename(dd) # get 'bar' from '/foo/bar/'

        date = dd.split('_')[0]
        date = date[:4]+'-'+date[4:6]+'-'+date[6:8] # Change to YYYY-MM-DD format.
        study_id = dd.split('_')[1]
        visit_id = dd.split('_')[2].split('/')[0]
        if replace_visit_prefix:
            # Change H02 -> V02; B00 -> V00, V11 -> V01 etc.
            visit_id = visit_id.replace(visit_id[0:2], "V0")

#        print dd, date, study_id, visit_id

        # Add row for Study ID if it does not exist.
        if (not study_id in report_df.index.tolist()):
            report_df = report_df.append(pandas.Series(name=study_id)) # index is Study ID
            report_df.fillna('', inplace=True)

        # Add column for visit ID if it does not exist.
        while True: # Change the visit ID until we find a visit ID tag that doesn't already contain a value.
            if (not visit_id in report_df.columns):
                report_df[visit_id] = ''
                report_df.fillna('', inplace=True)

            # Add data to cell.
            if (report_df.loc[study_id, visit_id] == ''):
                report_df.loc[study_id, visit_id] = date
                break
            else:
                if (visit_id[-1] != ")"):
                    visit_id += " (2)"
                else:
                    last_repeat_number = int(visit_id[-2])
                    visit_id = visit_id[:-3]+"("+str(last_repeat_number+1)+")"
            
    # Fill empty cells with ''.
    report_df.fillna('', inplace=True)

    # Sort the Visit IDs based on the last digit.
    try:
        cols_sorted = sorted(report_df.columns, key=lambda x:int(x[2]+x[1]))
    except ValueError:
        cols_sorted = sorted(report_df.columns, key=lambda x:int(x[2]))
    report_df = report_df[cols_sorted]

    # Sort the study IDs.
    report_df.sort_index(inplace=True)
    
    # Convert to dates
#    for col in report_df.columns:
#        report_df[col] = report_df[col].dt.date
    
    # Save report.
    report_df.to_csv(out_path)
    
    return

if (__name__ == '__main__'):
    # Create argument parser
    description = """Generate a report on the dates at which each patient had their V01, V02 scans etc., based on the directory names in /hpf/largeprojects/smiller/images/"""
    epilog = None
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("in_dir", help="Directory in which to search for and report visit dates.")
    parser.add_argument("out_path", help="Path to which report will be saved. Must have extension .csv.")    

    # Define keyword argument.
    parser.add_argument("-i", "--ignore_symlinks", help="Do not report on symbolic links. Still raises a warning.", action="store_true")
    parser.add_argument("-v", "--replace_visit_prefix", help="If different visit codes have different letters or first digits (e.g. for CND study, use B## for brain, H## for heart, and C## for combined, or V01, V11, V21 for repeated V#1 scans), treat all letters and first digits as the same.", action="store_true")
    
    # Print help if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Do stuff.
    visitIdReport(in_dir=args.in_dir, out_path=args.out_path, ignore_symlinks=args.ignore_symlinks, replace_visit_prefix=args.replace_visit_prefix)
