#!/usr/bin/env python

#### Python code to deidentify a DICOM file
## This script was written for use in a UNIX or UNIX-like operating system, which includes the UNIX 'file' function.
## This script requires the functions 'dcmdump' and 'dcmodify', part of the DICOM Toolkit (DCMTK). These functions must be avaiable in the PATH.
## For usage, do:
# $ python anonDicom.py -h



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

def anonDicom(in_path, name=None, level=2, modify_pid=False, print_only=False, force=False):
    # Check if file is a directory or is not DICOM
    if os.path.isdir(in_path):
        print "Error: Input path is a directory. To anonymize all files in a directory, use the '-r' flag."
        return

    # Check if file appears to be DICOM using the Unix command 'file'. This may fail for real DICOM files which are missing information in the header. In such cases, consider running anonDicom with the '-f --force' flag enabled, to skip this check.
    if (not force):
        cmd_file = 'file "%s"' % (in_path)
        output, errors = run_cmd(cmd_file, 0)
        if (not "DICOM medical imaging data" in output):
            err_msg = "Error: File appears to be non-DICOM or corrupted. Not anonymizing: %s" % (in_path)
            print err_msg
            return

    # Define list of tags to remove (excluding PatientName and PatientID)
    tags = {} # keys are tag types (e.g. LO, LT, PN); elems are lists of tags.
    
    type_list = ["AE", "AS", "CS", "DA", "DS", "LO", "LT", "PN", "SH", "SQ", "ST", "TM", "US"]

    for type in type_list:
        tags[type] = []

    if (level >= 1):
        tags["CS"].append("(0010,2298)") # Responsible Person Role
        tags["DA"].append("(0010,0030)") # Patient's Birth Date
        tags["LO"].append("(0010,0033)") # Patient's Birth Date in Alternative Calendar
        tags["LO"].append("(0010,0034)") # Patient's Death Date in Alternative Calendar
        tags["LO"].append("(0010,1000)") # Other Patient IDs
        tags["LO"].append("(0010,1040)") # Patient's Address
        tags["LO"].append("(0010,1050)") # Insurance Plan Identification (RET)
        tags["LO"].append("(0010,1090)") # Medical Record Locator
        tags["LO"].append("(0010,2299)") # Responsible Organization
        tags["LT"].append("(0008,4000)") # Identifying Comments (RET)
        tags["LT"].append("(0010,4000)") # Patient Comments
        tags["PN"].append("(0010,1001)") # Other Patient Names
        tags["PN"].append("(0010,1005)") # Patient's Birth Name
        tags["PN"].append("(0010,1060)") # Patient's Mother's Birth Name
        tags["PN"].append("(0010,2297)") # Responsible Person
        tags["SH"].append("(0010,2154)") # Patient's Telephone Numbers
        tags["SQ"].append("(0010,0050)") # Patient's Insurance Plan Code Sequence
        tags["SQ"].append("(0010,1002)") # Other Patient IDs Sequence
        tags["TM"].append("(0010,0032)") # Patient's Birth Time

    if (level >= 2):
        tags["AS"].append("(0010,1010)") # Patient's Age
        tags["CS"].append("(0010,0040)") # Patient's Sex
        tags["CS"].append("(0010,21A0)") # Smoking Status
        tags["DA"].append("(0010,21D0)") # Last Menstrual Date
        tags["DS"].append("(0010,1020)") # Patient's Size
        tags["DS"].append("(0010,1030)") # Patient's Weight
        tags["LO"].append("(0010,1080)") # Military Rank
        tags["LO"].append("(0010,1081)") # Branch of Service
        tags["LO"].append("(0010,2000)") # Medical Alerts
        tags["LO"].append("(0010,2110)") # Allergies
        tags["LO"].append("(0010,2150)") # Country of Residence
        tags["LO"].append("(0010,2152)") # Region of Residence
        tags["LO"].append("(0010,21F0)") # Patient's Religious Preference
        tags["LO"].append("(0032,1030)") # Reason for Study (RET)
        tags["LO"].append("(0038,0300)") # Current Patient Location
        tags["LO"].append("(0038,0400)") # Patient's Institution Residence
        tags["LO"].append("(0038,0500)") # Patient State
        tags["LT"].append("(0010,21B0)") # Additional Patient History
        tags["LT"].append("(0032,4000)") # Study Comments
        tags["PN"].append("(0008,0090)") # Referring Physician's Name
        tags["PN"].append("(0008,009C)") # Consulting Physician's Name
        tags["PN"].append("(0008,1048)") # Physician(s) of Record
        tags["PN"].append("(0008,1050)") # Performing Physician's Name
        tags["PN"].append("(0008,1060)") # Name of Physician(s) Reading Study
        tags["PN"].append("(0008,1070)") # Operators' Name
        tags["PN"].append("(0032,1032)") # Requesting Physician
        tags["PN"].append("(0040,0006)") # Scheduled Performing Physician's Name
        tags["PN"].append("(0040,2008)") # Order Entered By
        tags["PN"].append("(4008,010C)") # Interpretation Author (RET)
        tags["PN"].append("(4008,0114)") # Physician Approving Interpretation (RET)
        tags["SH"].append("(0008,0094)") # Referring Physician's Telephone Numbers
        tags["SH"].append("(0010,2160)") # Ethnic Group
        tags["SH"].append("(0010,2180)") # Occupation
        tags["SH"].append("(0040,2009)") # Order Enterer's Location
        tags["SH"].append("(0040,2010)") # Order Callback Phone Number
        tags["SQ"].append("(0008,0096)") # Referring Physician Identification Sequence
        tags["SQ"].append("(0008,009D)") # Consulting Physician Identification Sequence
        tags["SQ"].append("(0008,1049)") # Physician(s) of Record Identification Sequence
        tags["SQ"].append("(0008,1052)") # Performing Physician Identification Sequence
        tags["SQ"].append("(0008,1062)") # Physician(s) Reading Study Identification Sequence
        tags["SQ"].append("(0008,1072)") # Operator Identification Sequence
        tags["SQ"].append("(0010,0101)") # Patient's Primary Language Code Sequence
        tags["SQ"].append("(0010,0102)") # Patient's Primary Language Code Modifier Sequence
        tags["SQ"].append("(0010,1021)") # Patient's Size Code Sequence
        tags["SQ"].append("(0032,1031)") # Requesting Physician Identification Sequence
        tags["SQ"].append("(0040,000B)") # Scheduled Performing Physician Identification Sequence
        tags["ST"].append("(0008,0092)") # Referring Physician's Address
        tags["ST"].append("(0008,0116)") # Responsible Organization
        tags["US"].append("(0010,21C0)") # Pregnancy Status
        
    if (level == 3):
        tags["AE"].append("(0032,1021)") # Scheduled Study Location AE Title (RET)
        tags["AE"].append("(0040,0241)") # Performed Station AE Title
        tags["CS"].append("(0010,0022)") # Type of Patient ID
        tags["CS"].append("(0010,0035)") # Patient's Alternative Calendar
        tags["LO"].append("(0008,0080)") # Institution Name
        tags["LO"].append("(0008,1040)") # Institutional Department Name
        tags["LO"].append("(0010,0021)") # Issuer of Patient ID
        tags["LO"].append("(0032,1020)") # Scheduled Study Location (RET)
        tags["PN"].append("(4008,0119)") # Distribution Name (RET)
        tags["SH"].append("(0008,1010)") # Station Name
        tags["SH"].append("(0040,0010)") # Scheduled Station Name
        tags["SH"].append("(0040,0011)") # Scheduled Procedure Step Location
        tags["SH"].append("(0040,0242)") # Performed Station Name
        tags["SH"].append("(0040,0243)") # Performed Location
        tags["SQ"].append("(0008,0082)") # Institution Code Sequence
        tags["ST"].append("(0008,0081)") # Institution Address        
        
    # Define name tags to be overwritten.
    name_tags = {}
    if (name != None):
        name_tags = {}
        name_tags["PN"] = ["(0010,0010)"]
        if (modify_pid): # if opting to overwrite PatientID in addition to PatientName
            name_tags["LO"] = ["(0010,0020)"]

    # Define replacement tags
    a_text = "ANONYMIZED" # Replace text with this string
    a_age  = "000Y"
    a_date = "00000000"
    a_time = "000000.00"
    a_num  = "0"

    if (not print_only): # if anonymizing rather than reporting on PHI.
        ## Build dcmodify command.
    
        # Start command and set high-level options
        cmd_dcmodify = 'dcmodify -nb -ie'
    
        # Add command for modifying each tag, type by type.
        # AE
        for tag in tags["AE"]:
            cmd_dcmodify += ' -ma "%s=%s"' % (tag, a_text)
        # AS
        for tag in tags["AS"]:
            cmd_dcmodify += ' -ma "%s=%s"' % (tag, a_age)
        # CS
        for tag in tags["CS"]:
            cmd_dcmodify += ' -ma "%s=%s"' % (tag, a_text)
        # DA
        for tag in tags["DA"]:
            cmd_dcmodify += ' -ma "%s=%s"' % (tag, a_date)
        # DS
        for tag in tags["DS"]:
            cmd_dcmodify += ' -ma "%s=%s"' % (tag, a_num)
        # LO
        for tag in tags["LO"]:
            cmd_dcmodify += ' -ma "%s=%s"' % (tag, a_text)
        # LT
        for tag in tags["LT"]:
            cmd_dcmodify += ' -ma "%s=%s"' % (tag, a_text)
        # PN
        for tag in tags["PN"]:
            cmd_dcmodify += ' -ma "%s=%s"' % (tag, a_text)
        # SH
        for tag in tags["SH"]:
            cmd_dcmodify += ' -ma "%s=%s"' % (tag, a_text)
        # SQ # Only tag type that is deleted. All others are modified only.
        for tag in tags["SQ"]:
            cmd_dcmodify += ' -ea "%s"' % (tag)
        # ST
        for tag in tags["ST"]:
            cmd_dcmodify += ' -ma "%s=%s"' % (tag, a_text)
        # TM
        for tag in tags["TM"]:
            cmd_dcmodify += ' -ma "%s=%s"' % (tag, a_time)
        # US
        for tag in tags["US"]:
            cmd_dcmodify += ' -ma "%s=%s"' % (tag, a_num)
    
        ## Name tags
        if (name != None):
            for tag in name_tags["PN"]:
                cmd_dcmodify += ' -ma "%s=%s"' % (tag, name)
                
            if modify_pid:
                for tag in name_tags["LO"]:
                    cmd_dcmodify += ' -ma "%s=%s"' % (tag, name)
        else:
#            print 'WARNING: No name specified, PatientName will not be changed.'
            pass

        cmd_dcmodify += ' "%s"' % (in_path)
    
        # Execute anonymization command.
        output, errors = run_cmd(cmd_dcmodify, 0)
#        print errors

    else: # if merely reporting PHI. Do this in a separate script.
        # Write dcmdump command to find and print all PHI tags that exist in the header.
        tags_string = ''
        name_tags = {}
        name_tags["PN"] = ["(0010,0010)"]
        name_tags["LO"] = ["(0010,0020)"]
        for tag_list in [name_tags, tags]:
            for key, tags in tag_list.iteritems():
                for tag in tags:
                    tags_string += tag+'|'
        tags_string = tags_string.rstrip('|') # remove trailing '|'.
        cmd_dcmdump = 'dcmdump "%s" | grep -E "%s" | grep -vE "no value available|%s"' % (in_path, tags_string, a_text)

        # Execute dcmdump command.
        output, errors = run_cmd(cmd_dcmdump, 0)
#        print output
#        print errors
    return

if (__name__ == "__main__"):
    # Create argument parser
    description = """Overwrite identifying tags in a DICOM file with anonymous values."""
    epilog = """Notes:
- The PatientName tag can be overwritten with a user-specified subject ID by using the '-n' option. If no subject ID is specified, the PatientName tag will not be modified.
- Using the '-m' option, the PatientID tag will also be overwritten with the user-specified subject ID. If the '-n' and '-m' options are not both used, the PatientID tag will not be modified.

Examples:
# Anonymize a single DICOM file without modifying the PatientName or PateintID tags:
anonDicom.py 0123.dcm

# Anonymize a single DICOM file, changing the PatientName tag to "ABC123", leaving the PatientID tag unchanged:
anonDicom.py 0123.dcm -n ABC123 

# Anonmize a single DICOM file, changing the PatientName and PatientID tags to "ABC123":
anonDicom.py 0123.dcm -n ABC123 -m

# Check anonymization of a single DICOM file without modifying the file:
anonDicom.py 0123.dcm -p

# Anonyize all DICOM files in a directory:
anonDicom.py /path/to/my/directory -r 
 """
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("in_path", help="path to DICOM file to be anonymized, or directory containing DICOM files.")
    
    # Define optional arguments.
    parser.add_argument("-n", "--name", help="Subject ID to set PatientName tags to", type=str)
#    parser.add_argument("-b", "--backup", help="add an option for backup. Do not back up by default.")
    parser.add_argument("-l", "--level", type=int, default=2, choices=[1,2,3], help="Set degree of anonymization (1: directly identifying information such as patient name, birth date, address, phone number. 2 (default): indirectly identifying information such as weight, age, physicians etc. 3: information about institution which performed scan, such as address, department etc.)")
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
        print "Error: Input path does not exist:", args.in_path
        sys.exit()

    # Build list of files to anonymize.
    if (not args.recursive):
        file_path_list = [args.in_path]
    else:
        if (not os.path.isdir(args.in_path)):
            print "Error: Input path is not a directory, but the '-r' flag was used. Quitting."
            sys.exit()
        file_path_list = buildPathList(args.in_path)
    num_files = len(file_path_list)
    
    # Anonymize
    if (args.recursive):
        if (num_files == 0):
            print "Error: No files found in input directory."
        else:
            print "Attempting to anonymize "+str(len(file_path_list))+" files."
    for file_path in file_path_list:
        anonDicom(file_path, args.name, level=args.level, modify_pid=args.modify_pid, print_only=args.print_only, force=args.force)

# WHAT I SHOULD DO WITH THIS SCRIPT:
# Create a lookup table for all the DICOM tags which specifies the tag type, importance level.
# Another table should specifed the value to change each tag type to.
# I should then have two separate commands, one which performs the anonymization, and another
# which checks for and prints possible PHI. Both of these scripts will rely on the same table
# and replacement values, so the script which checks for PHI can be configured to only print
# the tag if the value stored in that tag is not the replacement value stored in the table.
