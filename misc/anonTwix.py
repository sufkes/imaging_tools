#!/usr/bin/env python

import os, sys, argparse

def getVal(lines, index):
    line=lines[index]
    if ('{' in line) and ('}' in line):
        val = line.split('{')[1].split('}')[0].strip()
        print val
    elif ('{' in lines[index+1]):
        val = lines[index+2].strip()
    else:
        raise ValueError('Cannot determine value of tag at line '+index)
    return val

def setVal(lines, index):
    return

def anonRDA(in_path, name=None, backup=False, verbose=False, print_only=False):
    """
    Parameters:
        in_path: path to TWIX file to modify
    Modifies header of TWIX file
    """

    ## Verify existence of file and open it.
    if (os.path.exists(in_path)):
        with open(in_path, 'rb') as fh:
            lines = fh.readlines()
    else:
        raise ValueError("Input path does not exist")

    ## Determine whether to anonymize or just print PHI-containing fields.
    if (print_only):
        # List of known PHI-containing fields; could be lots more.
#        tags_to_print = ['PatientName', 'PatientID', 'PatientSex', 'PatientBirthDate', 'PatientAge', 'PatientWeight']
#        for tag in tags_to_print:
#            for line in lines:
#                if tag in line:
#                    print line,
        pass
            
    else:
        ## Create backup.
        if (backup):
            backup_path = in_path+'.bak'
            while (os.path.exists(backup_path)):
                backup_path += ".bak"
            with open(backup_path, 'w') as fh:
                fh.writelines(lines)
    
        ## Modify header tags.
        # Build list of header tags to modify.
        
#        tags = {'PatientSex':'A', 'PatientBirthDate':'00000000', 'PatientAge':'999Y', 'PatientWeight':'0.000000'}
#        if (name != None):
#            tags['PatientName'] = name
#            tags['PatientID'] = name
#        else:
            # Find PatientName.
            
        for index in range(len(lines)):
            line = lines[index]
            if ('ParamString."tPatientName"' in line) or ('ParamString."PatientName"' in line):
                print getVal(lines, index)
                #pass
#            for line in lines:
#                if (line[:11] == 'PatientName'):
#                    PatientName = ": ".join(line.split(": ")[1:]).split("\r\n")[0]
#                    tags['PatientID'] = PatientName
#                    break
    
        # Modify header lines. Only modify fist line matching header name.
#        for tag, tag_val_new in tags.iteritems():
#            for line_num in range(len(lines)):
#                line = lines[line_num]
#                if tag in line:
#                    tag_val_old = ": ".join(line.split(": ")[1:]).split("\r\n")[0]
#                    lines[line_num] = line.replace(tag_val_old, tag_val_new) 
    
        # Save modified file.
        with open(in_path+'new', 'wb') as fh:
            fh.writelines(lines)

    return

if (__name__ == "__main__"):
    # Create argument parser.
    description = """
    Anonymize TWIX file header.
    """
    parser = argparse.ArgumentParser(description=description)
    
    # Define positional arguments.
    parser.add_argument("in_path", help="Path to Twix file to be anonymized", type=str)
        
    # Define optional arguments.
    parser.add_argument("-n", "--name", type=str, action="store", help="Name to change PatientName and PatientID to. By default, PatientID is set to the current PatientName, and PatientName is unchanged.")
    parser.add_argument("-b", "--backup", action="store_true", help="Create backup of the file.")
    parser.add_argument("-v", "--verbose", action="store_true", help="[NONFUNCTIONAL] Print the changes made to the header file.")
    parser.add_argument('-p', '--print-only', help='Print PHI-containing tags. Do not anonymize.', action='store_true')
    
    # Print help message if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()
    
    # Parse arguments.
    args = parser.parse_args()
    
    # Perform check
    anonRDA(args.in_path, name=args.name, backup=args.backup, verbose=args.verbose, print_only=args.print_only)
