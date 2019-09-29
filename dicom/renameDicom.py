#!/usr/bin/python
#  Python code to rename dicom files sent from any scanner
#  Adapted from renameDicom.pl
#  Created to allow for great portability and flexibility between scanners
#       Checks dicom header if Philips scanner
#  Increased options with respect to output extension and anonymization 

#    File Name:  process_fmri.py
#
# OUTPUT STRUCTURE/home/researchPACS/dcmsrvr/datadump/MR.1.2.840.113619.2.312.4120.8419247.14553.1343646753.9
#	base_dir / series_dir / file_name
#
#	base_dir = <StudyDate>_<SubjectID>
#	series_dir = <SeriesNumber>_<ProtocolName>
#	file_name = I<InstanceNumber>.dcm
#		- InstanceNumber is padded to 4 digits
#
# LAST REVISION 
#	12/05/01 - WL - initial creation 
#   12/05/11 - WL - Fixed modify when moving, added valid char check
#   12/05/15 - WL - Added echo number (0018,0086) check in case of duplicate file
#                   Restore clobber functionality
#   12/05/30 - WL - Added error checking for echo number field 
#                   Modified to match naming convention used by Josh
#   12/08/07 - WL - Removed StudyDescription field, it's unused anyway

import os
import string
import shlex, subprocess
import datetime
from optparse import OptionParser, Option, OptionValueError


program_name = 'renameDicom.py'

# Defining dcm header names and tags
lut_tag  = {}
lut_tag['StudyDate'] = '0008,0020'
lut_tag['StudyTime'] = '0008,0030'
lut_tag['SeriesNum'] = '0020,0011'
#lut_tag['StudyDescription'] = '0008,1030'
# lut_tag['ProtocolName'] = '0018,1030'
lut_tag['InstanceNumber'] = '0020,0013'
lut_tag['SeriesDescription'] = '0008,103e'
lut_tag['PatientsName'] = '0010,0010'
lut_tag['StudyInstanceUID'] = '0020,000d'
# SKU 2019-01-11: Include SeriesInstanceUID in series folder name to avoid putting different series in same folder.
lut_tag['SeriesInstanceUID'] = '0020,000e'

# What characters are valid in a file name?
# Anything not in this group is removed
valid_chars = '-_%s%s' % (string.ascii_letters, string.digits)

lut_value = {}     # Blank structure for storing header values

def run_cmd(sys_cmd, debug, verbose):
# one line call to output system command and control debug state
    if verbose:
        print sys_cmd
    if not debug:
        p = subprocess.Popen(sys_cmd, stdout = subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = p.communicate()
        return output, errors
    else:
        return '','' 
        
if (__name__ == '__main__'):
    usage = "Usage: "+program_name+" <options> dir_input fname_dcm dir_output\n"+\
            "   or  "+program_name+" -help\n" +\
            " For directories use: \n" +\
            " for file in <dir>; do renameDicom.py . ${file} <dir_output>; done"
    parser = OptionParser(usage)
    parser.add_option("-c","--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    parser.add_option("-a","--anon", action="store_true", dest="anon",
                        default=0, help="make dicom anonymous")
    parser.add_option("--cfg", type="string", dest="fname_cfg",
                        default="DcmServer.cfg", help="List of fields to wipe [default = DcmServer.cfg]")
    parser.add_option("-m","--move", action="store_true", dest="move",
                        default=0, help="move instead of copy")
    parser.add_option("-v","--verbose", action="store_true", dest="verbose",
                        default=0, help="Verbose output")
    parser.add_option("-d","--debug", action="store_true", dest="debug",
                        default=0, help="Run in debug mode")
    parser.add_option("-e","--extension", type="string", dest="extension",
                        default="dcm", help="File extension [default = dcm]")
    parser.add_option("-l","--logile", type="string", dest="logfile",
                        default="", help="Log file location [default = none]")

# # Parse input arguments and store them
    options, args = parser.parse_args()     
    if len(args) != 3:
        parser.error("incorrect number of arguments")
    dir_input, fname_dcm, dir_output = args

    if options.move:   # Determine what type of operation to carry out
        operation = 'mv'
    else:
        operation = 'cp'
    
    # Start line to add to the log.
    time_stamp = str(datetime.datetime.now()).split('.')[0]  # grap timestamp, but drop ms; moved to start of script.
    line_log = '' # Add log lines to this string; add them to the log in one step.

    # Check that file exists
    if not os.path.exists('%s/%s' % (dir_input, fname_dcm)):
        raise SystemExit, 'DCM file does not exist - %s/%s' % \
            (dir_input, fname_dcm)
    
    # Check scanner type (0008,0070 - Manufacturer tag)
    cmd_dcmdump = 'dcmdump %s/%s | grep 0008,0070' % (dir_input, fname_dcm)
    output, errors = run_cmd(cmd_dcmdump,0, 0)
    if string.lower(output).find('philips') > 0:
        scanner_type = 'philips'
    elif string.lower(output).find('siemens'):
        scanner_type = 'siemens'
    elif string.lower(output).find('ge'):
        scanner_type = 'ge'
    else:
        scanner_type = 'other'

    for tag_name in lut_tag:
        # Dump dcmheader and grap header line
        cmd_dcmdump = 'dcmdump %s/%s | grep %s' % (dir_input, fname_dcm, lut_tag[tag_name])
        output, errors = run_cmd(cmd_dcmdump, 0, 0)

        # Check if Philips type dicom and looking for InstanceNumber
        if (scanner_type == 'philips') and (tag_name == 'InstanceNumber'):
            valid_line_found = False
            for line in output.split('\n'): # split output into separate lines
                if (line != '') and (line.find('[0]') == -1): # Ignore all lines without 'real' instance information
                    output = line # *use last line that has a nonzero value. Is this desired behaviour?
                    valid_line_found = True

        # 2019-02-08 Steven Ufkes - general check for multiple tags
        elif (len(output.split('\n')) > 2):
            line_log += time_stamp + " - WARNING: Multiple values found for tag '"+tag_name+"' in file "+dir_input+"/"+fname_dcm+'\n'
            valid_line_found = False
            for line in output.split('\n'):
                if (not line=='') and (line.find('(no value available)') == -1): # Ignore all lines without 'real' instance information; Use first line found with data.
                    output = line
                    valid_line_found = True
                    break

        # 2019-02-08 Steven Ufkes - Assume exactly one line with the tag will be found. 
        # If the line is missing, this error will be caught when the tag_value is 
        # extracted from the line.
        else:
            valid_line_found = True

        # Define error message for log it tag value can't be found.
        missing_tag_msg = time_stamp + " - WARNING: Missing DICOM tag '"+tag_name+"' in file "+dir_input+"/"+fname_dcm+'\n'
        # Define tag value to use if it cannot be determined.
        backup_tag_value = "Unknown"+tag_name

        if valid_line_found:
            try:
                tag_value = output.split('[')[1].split(']')[0]    # Grab value between [ ]
            except IndexError: # if line with correct tag number and containing '[something]' is not found.
                line_log += missing_tag_msg
                tag_value = backup_tag_value
        else:
            line_log += missing_tag_msg
            tag_value = backup_tag_value

        # Initial replacing of bad characters into useful partitions
        tag_value = tag_value.replace(' ','-')   # replace spaces with -  12/05/30 - WL - Custom
        tag_value = tag_value.replace('.','-')   # replace . with -       12/05/30 - WL - Custom
        tag_value = tag_value.replace('/','-')   # replace / with -
        tag_value = tag_value.replace('\'','-')   # replace \ with - # SKU 2019-02-01: I think this replaces ' with -
        tag_value = tag_value.replace('\\','_')  # SKU 2019-02-08 - replace \ with -
        tag_value = tag_value.replace('*','s')   # replace * with s
        tag_value = tag_value.replace('?','q')   # replace ? with q
        tag_value = ''.join(c for c in tag_value if c in valid_chars) # scrub bad characters
        lut_value[tag_name] = tag_value
#        print tag_name + ' = ' , tag_value
    
    # Creating output directories and namesq
		# 12/05/30 - WL - Custom
    dir_base = '%s_%s_%s' % \
        (lut_value['StudyInstanceUID'], lut_value['StudyDate'], lut_value['PatientsName']) 

    dir_base_full = '%s/%s' % (dir_output,dir_base)

    if not os.path.exists(dir_base_full):
        cmd_mkdir = 'mkdir %s/%s' % (dir_output, dir_base)
        output, errors = run_cmd(cmd_mkdir, options.debug, options.verbose)

    dir_series = '%03.d-%s-UID_%s' % (int(lut_value['SeriesNum']), lut_value['SeriesDescription'], lut_value['SeriesInstanceUID'])
    dir_series_full = '%s\%s' % (dir_base_full, dir_series)
    if not os.path.exists(dir_series_full):
        cmd_mkdir = 'mkdir %s/%s/%s' % (dir_output, dir_base, dir_series)
        output, errors = run_cmd(cmd_mkdir, options.debug, options.verbose)

        
    # iterative process to determine file name
    # if lut_value['ProtocolName'] == 'Phoenix_Document':
    #     fname_out = 'Series_%s' % (int(lut_value['InstanceNumber']),)
    #     extension = 'SR'
    if False:
    	pass
    else:
	# 12/05/30 - WL - Custom 
#        fname_out = '%04.d' % (int(lut_value['InstanceNumber']))

        # SKU 2019-01-11: Try simply not renaming the dicom files (except appending A's in case of duplicates).
        fname_out = fname_dcm

        extension = options.extension
        
    full_out = '%s/%s/%s/%s.%s' % (dir_output, dir_base, dir_series, fname_out, extension)
    
# If Siemens scanner, check for possibility of Mag/Ph type output, where files will have the same 
# instance number (leading to the same destination file name) but different echo numbers (0018,0086)
# This problem does not appear to exist on GE and Philips scanners (which have different echo #s and
# instance #s
#
# SKU 2019-01-11: This section seems sketchy, and is likely the cause of the problem of files being overwritten.
# I have commented out the section below. Now, images should all be named in the same way regardless of
# echo numbers. There should not be anyway for files to be removed or overwritten. In case of duplicate 
# destination files, A's should be appended to the most recently added file in every case.
# 
#    if scanner_type  == 'siemens':
#        # get source echo number
#        cmd_dcmdump = 'dcmdump %s/%s | grep 0018,0086' % (dir_input, fname_dcm)
#        output, errors = run_cmd(cmd_dcmdump, 0, 0)
#        if output != '':
#            echo_number_source = int(output.split('[')[1].split(']')[0])
#        else:
#            echo_number_source = 1
#        if echo_number_source > 1:   # append echo number if it's >1
#            # if non-suffixed file exists (ie. echo number = 1) then add suffix
#            if os.path.exists(full_out):
#                if options.verbose:
#                    print 'Renaming destination file - first echo'
#                fname_new = '%s_01.%s' % (os.path.splitext(full_out)[0], extension)
#                cmd_mvdest = 'mv %s %s' % (full_out, fname_new) # SKU: This could overwrite files.
#                output, errors = run_cmd(cmd_mvdest, options.debug, options.verbose)
#                if errors == '':
#                    cmd_rmorig = 'rm %s' % (full_out) # SKU: Can't see how 'full_out' could ever exist at this point
#                    output, errors = run_cmd(cmd_rmorig, options.debug, options.verbose)
#            fname_out = '%s_%02.d' % (fname_out, echo_number_source)
#            full_out = '%s/%s/%s/%s.%s' % \
#                (dir_output, dir_base, dir_series, fname_out, extension)
#        else:
#            # if new file
#            if not os.path.exists(full_out):
#                dir_full_out = '%s/%s/%s' % (dir_output, dir_base, dir_series)
#                list_existing = str(os.listdir(dir_full_out))
#                # check if basic filename already exists
#                if list_existing.find(fname_out)>-1: 
#                    fname_out = '%s_01' % (fname_out,)
#                full_out = '%s/%s/%s/%s.%s' % \
#                    (dir_output, dir_base, dir_series, fname_out, extension)
 
## 2019-01-25 - Steven Ufkes - This section is insidiuosly bugged. It will successfully find a 
## filepath which does not exist. However, it is possible that the file path will exist at the 
## time when the mv command is executed. The Unix mv command was being called in clobber mode,
## and would occasionally overwrite a file that was created between the time at which the
## destination filepath was chosen and the time at which mv was called
#    if not options.clobber:    # if not overwriting then must find unique filename 
#        while os.path.exists(full_out):
#            fname_out = fname_out + 'A'
#            full_out = '%s/%s/%s/%s.%s' % (dir_output, dir_base, dir_series, fname_out, extension)
    
#    time_stamp = str(datetime.datetime.now()).split('.')[0]  # grap timestamp, but drop ms; moved to start of script.

## 2019-01-25 - Steven Ufkes - This section was rejigged per the comments above.
#    cmd_mvdcm = '%s %s/%s %s' % \
#        (operation, dir_input, fname_dcm, full_out)
#    output, errors = run_cmd(cmd_mvdcm, options.debug, options.verbose)     # Copying DCM with new name
#### 2019-01-25 - Steven Ufkes (SKU) - Start of newly added block
    if options.clobber:
        cmd_mvdcm = '%s %s/%s %s' % \
            (operation, dir_input, fname_dcm, full_out)
        overwrite_completed = os.path.exists(full_out) # *Could fail to detect overwrite due to race condition.
        output, errors = run_cmd(cmd_mvdcm, options.debug, options.verbose)     # Copying DCM with new name
    else:
        # To avoid the bug described above (2019-01-25), write a temporary file
        # to the destination path if a file doesn't exist. This will prevent 
        # multiple processes from writing to the same destination at the same time,
        # because only one process will be able to write the temporary file if the 
        # os.O_EXCL flag is used.
        operation_completed = False
        cmd_overwrite_attempt = None
        while (not operation_completed):
            try:
                cmd_mvdcm = '%s %s/%s %s' % (operation, dir_input, fname_dcm, full_out) # overwrite temp file with dicom file
		os.open(full_out, os.O_CREAT|os.O_EXCL) # create empty temp file at destination path
                output, errors = run_cmd(cmd_mvdcm, options.debug, options.verbose) # Copying DCM with new name
                operation_completed = True
            except OSError:
                if (cmd_overwrite_attempt == None): # If this is the first overwrite attempt for this file.
                    cmd_overwrite_attempt = cmd_mvdcm # Store the first overwrite attempt command for the log.
                fname_out = fname_out + 'A'
                full_out = '%s/%s/%s/%s.%s' % (dir_output, dir_base, dir_series, fname_out, extension)

#### 2019-01-25 - Steven Ufkes (SKU) - End of newly added block
    line_log += time_stamp + ' - ' + cmd_mvdcm + '\n'
    if (not options.clobber) and (cmd_overwrite_attempt != None): # report overwrite attempt in log (attempt to move file to existing path).
        line_log += time_stamp + ' - ' + "WARNING: overwrite attempted: '" + cmd_overwrite_attempt +"'\n"
    elif (options.clobber) and (overwrite_completed):
        line_log += time_stamp + ' - ' + "WARNING: overwrite completed: '" + cmd_mvdcm +"'\n"

    # Write executed commands and warnings to the log file.
    if (not options.debug):
        if (not options.logfile==''):   # log file
            file_log = open(options.logfile,'a')
            file_log.write(line_log)        
