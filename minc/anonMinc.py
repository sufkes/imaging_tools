#!/usr/bin/env python

#### Python code to deidentify a MINC file
## This script was written for use in a UNIX or UNIX-like operating system, which includes the UNIX 'grep' function.
## This script requires the functions 'mincheader' and 'minc_modify_header', part of the MINC Toolkit (DCMTK). These functions must be avaiable in the PATH.
## For usage, do:

import os, sys, argparse
import subprocess

def run_cmd(sys_cmd, verbose):
    if verbose:
        print sys_cmd
    p = subprocess.Popen(sys_cmd, stdout = subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = p.communicate()
    return output, errors

def buildPathList(in_dir): # get paths of all files in in_dir
    pathList = []
    for dirpath, dirnames, filenames in os.walk(in_dir):
        for filename in filenames:
            pathList.append(os.path.join(dirpath, filename))

    # Sort path list
    pathList = sorted(pathList, key=lambda p: (os.path.sep not in p, p))
    return pathList

def anonMinc(in_path, name=None, modify_pid=False, print_only=False, force=False):
    # Check if file is a directory or is not DICOM
    if os.path.isdir(in_path):
        raise ValueError('Input path is a directory')

    # Check if file appears to be MINC using the Unix command 'file'. This may fail for real MINC files which are missing information in the header. In such cases, consider running anonDicom with the '-f --force' flag enabled, to skip this check.
    # This test will probably result in false positives for non-minc files.
    if (not force):
        cmd_file = 'file "%s"' % (in_path)
        output, errors = run_cmd(cmd_file, 0)
        if (not "Hierarchical Data Format" in output) and (not "NetCDF Data Format data" in output): # These are the two responses I've seen from the 'file' for MINC files.
            print "Warning: File appears to be non-DICOM or corrupted. Not anonymizing: "+in_path
            return

    # Define list of tags to remove (excluding PatientName and PatientID)
    name_tags=['patient:full_name',
               'patient:identification']
    
    tags=['study:admitting_diagnosis',
          'study:attending_physician',
          'study:operator',
          'study:performing_physician',
          'study:radiologist',
          'study:referring_physician',
          'study:study_id',
          'patient:address',
          'patient:age',
          'patient:birthdate',
          'patient:insurance_id',
          'patient:other_ids',
          'patient:other_names',
          'patient:sex',
          'patient:size',
          'patient:weight']
    
    if (not print_only): # if anonymizing rather than reporting on PHI.
        ## Figure out which tags are actually used in the MINC header (since we only want ot insert anonymous values into tags that actually exist).
        # Build mincheader command.
        cmd_mincheader = 'mincheader '+in_path+" | grep -E 'patient\:|study\:'"

        # Execute mincheader command.
        output_mincheader, errors_mincheader = run_cmd(cmd_mincheader, 0)

        tags_used = []
        for tag in tags:
            if (tag in output_mincheader):
                tags_used.append(tag)
                
        ## Build minc_modify_header command.    
        # Start command and set high-level options
        cmd_mincmodify = "minc_modify_header "+in_path
    
        # Add command for modifying each tag, type by type.
        for tag in tags_used:
            cmd_mincmodify += " -sinsert '"+tag+"='"

        # Modify full_name and identification tags.
        if (name != None): # if a new name was specified.
            cmd_mincmodify += " -sinsert 'patient:full_name="+name+"'"
            if modify_pid:
                cmd_mincmodify += " -sinsert 'patient:identification="+name+"'"
        
        # Execute anonymization command.
        output_mincmodify, errors_mincmodify = run_cmd(cmd_mincmodify, 0)
        if (len(output_mincmodify.strip()) > 0):
            print output_mincmodify,
        if (len(errors_mincmodify.strip()) > 0):
            print errors_mincmodify,

    else: # if merely reporting PHI. Do this in a separate script.
        # Write mincheader command to find and print all PHI tags that exist in the header.
        cmd_mincheader = 'mincheader '+in_path+' | grep -E '
        cmd_mincheader += "'"+tags[0]
        for tag in tags[1:] + name_tags:
            cmd_mincheader += '|'+tag.replace(':', '\:')
        cmd_mincheader += "'"

        # Execute mincheader command.
        output_mincheader, errors_mincheader = run_cmd(cmd_mincheader, 0)
        if (len(output_mincheader.strip()) > 0):
            print output_mincheader,
        if (len(errors_mincheader.strip()) > 0):
            print errors_mincheader,
    return

if (__name__ == "__main__"):
    # Create argument parser
    description = """Delete information in identifying tags in a MINC file with anonymous values."""
    epilog = """This script has hardly been tested at all. Use with caution.
 """
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("in_path", help="path to MINC file to be anonymized, or directory containing DICOM files.")
    
    # Define optional arguments.
    parser.add_argument("-n", "--name", help="Subject ID to set PatientName tags to", type=str)
#    parser.add_argument("-b", "--backup", help="add an option for backup. Do not back up by default.")
    parser.add_argument('-m', "--modify_pid", help="Change PatientID to specified Subject ID. Default: False", action="store_true")
    parser.add_argument('-p', '--print-only', help='Print PHI-containing tags. Do not anonymize.', action='store_true')
    parser.add_argument('-r', '--recursive', action='store_true', help='if in_path is a directory, find and anononymize all files in that directory.')
    parser.add_argument('-f', '--force', action='store_true', help="If the Unix 'file' command does not identify the input file as DICOM, still attempt to anonymize the file. By default, files not identified as DICOM are skipped.")

    # Print help if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Check that input path exists
    if (not os.path.exists(args.in_path)):
        raise ValueError("Input path does not exist:", args.in_path)

    # Build list of files to anonymize.
    if (not args.recursive):
        file_path_list = [args.in_path]
    else:
        if (not os.path.isdir(args.in_path)):
            raise ValueError("Input path is not a directory, but the '-r' flag was used.")
        file_path_list = buildPathList(args.in_path)
    num_files = len(file_path_list)
    
    # Anonymize
    if (args.recursive):
        if (num_files == 0):
            raise ValueError("No files found in input directory.")
        else:
            print "Attempting to anonymize "+str(len(file_path_list))+" files."

    # Anonymize
    for file_path in file_path_list:
        anonMinc(file_path, name=args.name, modify_pid=args.modify_pid, print_only=args.print_only, force=args.force)

