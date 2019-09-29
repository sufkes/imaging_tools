#!/usr/bin/env python

import os, sys
import re
import argparse

# Load color printing function.
sufkes_git_repo_dir = "/Users/steven ufkes/scripts" # change this to the path to which the sufkes Git repository was cloned.
sys.path.append(os.path.join(sufkes_git_repo_dir, "misc"))
from Color import Color

def makeCodeBreaker(code_breaking_file_path):
    """
    Parameters:
        code_breaking_file_path: path to code-breaking CSV file in format [MRN, SubjectID, other irrelevant information]
    Returns:
        code_breaker: dict with MRNs as keys and SubjectIDs as values
    """
    code_breaker = {}
    with open(code_breaking_file_path, 'rb') as handle:
        for line in handle.readlines():
            columns = line.split(',')
            mrn = columns[0]
            id = columns[1]
            if (not mrn == 'NA'):
                if (mrn in code_breaker.keys()) and (id != code_breaker[mrn]):
                    raise ValueError("MRN corresponds to multiple subject IDs: "+mrn)
                code_breaker[mrn] = id
    return code_breaker

def findMRN(lines):
    """
    Parameters:
        lines: list, lines in MRI report
    Returns: 
        mrn: str, medical record number found in file
    """
    mrn = None
    for line in lines:
        if ('Medical Record Number:' in line):
            candidate_mrn = line.split('Medical Record Number:')[1].split('Birth Date:')[0].strip()
            if (mrn == None):
                mrn = candidate_mrn
            elif (mrn != candidate_mrn):
                raise ValueError('Mulitple unique MRNs found in MRI report')
    if (mrn == None):
        raise ValueError('No MRN found in MRI report.')
    return mrn

def anonLines(mri_report_lines, id):
    anon_lines = [] # fill with anonymized lines.
    num_sex_replace = 0

    # Create map of PHI keywords, the text which immediately follows them, their anonymous replacement value, and store the number of times the keyword is changed.
    default_replace_word = '[ANONYMIZED]'
    phi_keywords = {'Patient Name:':
                        {'next_word':'Sex:', 
                         'replace_word':id, 
                         'num_found':0}, 
                    'Sex:':
                        {'next_word':'Medical Record Number:', 
                         'replace_word':default_replace_word, 
                         'num_found':0}, 
                    'Medical Record Number:':
                        {'next_word':'Birth Date:', 
                         'replace_word':default_replace_word, 
                         'num_found':0}, 
                    'Birth Date:':
                        {'next_word':'Accession Number:', 
                         'replace_word':default_replace_word, 
                         'num_found':0}}


    sex_keywords = ['girl', 'boy', 'male', 'female']
#    unhandled_phi_keywords = ['woman', 'man', 'lady', 'she', 'her', 'daughter', 'sister', 'man', 'he', 'his', 'him', 'son', 'brother']

    for line in mri_report_lines:
        anon_line = line

        # Replace values for PHI keywords in standard locations.
        for keyword, keyword_info in phi_keywords.iteritems():
            next_word = keyword_info['next_word']
            replace_word = keyword_info['replace_word']
            if (keyword in anon_line):
                # Doesn't handle case in which keyword appears multiple times in the same line
                if (anon_line.count(keyword) > 1):
                    raise ValueError("PHI Keyword '"+keyword+"' appears multiple times in line. This script cannot handle this case.")

                anon_line = next_word.join([anon_line.split(keyword)[0] + keyword + ' '+replace_word+' '] + anon_line.split(next_word)[1:])
                phi_keywords[keyword]['num_found'] += 1
        
        # Remove indications of sex in free-text areas. 
        for sex_keyword in sex_keywords:
            if sex_keyword in anon_line.lower():
                anon_line = re.sub(sex_keyword, "[ANONYMIZED]", anon_line, flags=re.I)
                num_sex_replace += 1

#        # Report unhandled words which might be PHI:
#        for keyword in unhandled_phi_keywords:
#            if keyword in line:
#                print Color.red+"Warning:"+Color.end+" Did not remove word '"+Color.green+keyword+Color.end+"' from line: "
#                print re.sub(keyword, Color.green+keyword+Color.end, anon_line, flags=re.I)

        # Replace stripped newline characters:
        if (line[-1] == '\n') and (anon_line[-1] != '\n'):
            anon_line += '\n'

        # Add anonymized line to list of anonymized lines
        anon_lines.append(anon_line)
        

    for keyword, keyword_info in phi_keywords.iteritems():
        num_found = keyword_info['num_found']
        print str(num_found) + " anonymizations of keyword '"+keyword+"'"
    print str(num_sex_replace) + " anonymizations of sex keywords : "+str(sex_keywords)
    return anon_lines

def reportChanges(phi_lines, anon_lines):
    unhandled_phi_keywords = ['woman', 'man', 'lady', 'he', 'she', 'her', 'daughter', 'sister', 'man', 'his', 'him', 'son', 'brother']

    # Report unhandled keywords
    
#    for anon_line in anon_lines:
    anon_all = '\n'.join(anon_lines)
    keywords_missed = set()
    colored_line = anon_all
    for keyword in unhandled_phi_keywords:
        if keyword in anon_all:
            keywords_missed.add(keyword)
            colored_line = re.sub(keyword, Color.green+keyword+Color.end, colored_line, flags=re.I)
    if (len(keywords_missed) > 0):
        print Color.red+"Warning:"+Color.end+" Did not remove words "+str(list(keywords_missed))+" from line:"
        print colored_line

    # Report all changed lines
    for line_num in range(len(phi_lines)):
        phi_line = phi_lines[line_num]
        anon_line = anon_lines[line_num]
        if (phi_line != anon_line):
            print ########## Change ##########
            print "< " + phi_line
            print "---"
            print "> " + anon_line
            
    return

def anonMRIReport(mri_report_path, out_path, code_breaking_file_path, force=False):
    # Read MRI report and convert it to a list of strings.
    with open(mri_report_path, 'rb') as handle:
        mri_report_lines = handle.readlines()

    # Create a dict which maps MRN to SubjectID.
    code_breaker = makeCodeBreaker(code_breaking_file_path)

    # Find the MRN for the MRI report.
    mrn = findMRN(mri_report_lines)

    # Get the SubjectID for the MRI report.
    if (not mrn in code_breaker):
        raise ValueError('MRN in MRI report is not in the code breaking file: '+mrn)
    id = code_breaker[mrn]
    
    # Remove PHI from the MRI 
    anon_lines = anonLines(mri_report_lines, id)

    # Report changes
    reportChanges(mri_report_lines, anon_lines)

    # Write anonymized MRI report
    if (not force) and os.path.exists(out_path):
        raise ValueError('Output path exists')
    else:
        with open(out_path, 'wb') as handle:
            handle.writelines(anon_lines)

    return

if (__name__ == '__main__'):
    # Ad hoc script for anonymizing CND MRI reports for Steph.
    # Assume a given MRI report file is for exactly one patient.

    # Handle arguments
    mri_report_path = '/Users/steven ufkes/Documents/miller/image_processing/CND/mri_report_example_delete.txt'

    description = """Ad hoc script for anonymizing CND MRI reports for Steph.
    Assume a given MRI report file is for exactly one patient."""

    parser = argparse.ArgumentParser(description=description)
    
    # Define positional arguments.
    parser.add_argument("source_path", help="path to MRI report to be anonymized", type=str)
    parser.add_argument("dest_path", help="path to save anonymized MRI report to", type=str)
        
    # Define optional arguments.
    parser.add_argument('-c', '--code_path', help='path to code-breaking CSV file', type=str, default='/Users/steven ufkes/Documents/miller/image_processing/CND/cnd_key_delete.csv')
    parser.add_argument("-f", "--force", action="store_true", help="overwrite destination path if it exists")
#    parser.add_argument("-v", "--verbose", action="store_true", help="For each file in path_1, print whether file is exactly duplicated, only PixelData duplicated, or not duplicated. By default, the files in path_1 which have neither an exact duplicate nor a PixelData duplicate are printed.")
    
    # Print help message if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()
    
    # Parse arguments.
    args = parser.parse_args()

    # Anonymize MRI report
    anonMRIReport(args.source_path, args.dest_path, args.code_path, force=args.force)
