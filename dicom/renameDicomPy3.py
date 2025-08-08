#!/usr/bin/env python
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
#   25/01/22 - ML - Converted the entire script to run in Python3. Added logging facilities. Added CLI interface and options.
import logging
import re
import shutil
import pydicom
import subprocess

from argparse import ArgumentParser
from pathlib import Path

logger = logging.getLogger(__name__)

# Defining dcm header names and tags
# SKU 2019-01-11: Include SeriesInstanceUID in series folder name to avoid putting different series in same folder.
lut_tag = {'StudyDate': '0008,0020', 'StudyTime': '0008,0030', 'SeriesNumber': '0020,0011',
           'StudyDescription': '0008,1030', 'ProtocolName': '0018,1030', 'InstanceNumber': '0020,0013',
           'SeriesDescription': '0008,103e', 'PatientName': '0010,0010', 'StudyInstanceUID': '0020,000d',
           'SeriesInstanceUID': '0020,000e'}
lut_value = {}  # Blank structure for storing header values


def run_cmd(sys_cmd, debug, verbose):
    # one line call to output system command and control debug state
    if verbose:
        logger.debug(sys_cmd)
    if not debug:
        p = subprocess.Popen(sys_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
        stdout, stderr = p.communicate()
        return stdout, stderr
    else:
        return '', ''


if __name__ == '__main__':
    parser = ArgumentParser(
        description="For directories use: \n for file in <dir>; do renameDicom.py . ${file} <dir_output>; done")
    parser.add_argument("dir_input", type=Path, help="input directory")
    parser.add_argument("fname_dcm", type=Path, help="input file name")
    parser.add_argument("dir_output", type=Path, help="output directory")
    parser.add_argument("-c", "--clobber", action="store_true", dest="clobber",
                        help="overwrite output file")
    parser.add_argument("--cfg", type=Path, dest="fname_cfg", default="DcmServer.cfg",
                        help="List of fields to wipe [default = DcmServer.cfg]")
    parser.add_argument("-m", "--move", action="store_true", dest="move",
                        help="move instead of copy")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                        help="Verbose output")
    parser.add_argument("-d", "--debug", action="store_true", dest="debug",
                        help="Run in dry run mode")
    parser.add_argument("-e", "--extension", type=str, dest="extension", default="dcm",
                        help="File extension [default = dcm]")
    parser.add_argument("-l", "--logile", type=Path, dest="logfile", default=Path("/dev/null"),
                        help="Log file location [default = none]")

    # Parse input arguments and store them
    options = parser.parse_args()

    # Logging setup
    FORMAT = '%(asctime)s - %(lineno)d - %(levelname)s - %(message)s'
    if options.verbose:
        LEVEL = logging.DEBUG
    else:
        LEVEL = logging.INFO
    if options.logfile == Path("/dev/null") or not options.logfile.parent.is_dir():
        logging.basicConfig(level=LEVEL, format=FORMAT)  # Stream if log is not possible
    else:
        logging.basicConfig(filename=options.logfile.resolve(), level=LEVEL, format=FORMAT)

    dcm_path = options.dir_input.resolve() / options.fname_dcm

    if not dcm_path.is_file():
        raise ValueError(f'DCM file does not exist - {str(dcm_path)}')
    logger.debug(f'DCM file found - {str(dcm_path)}')
    dcm_data = pydicom.dcmread(dcm_path, stop_before_pixels=True)

    # Check scanner type (0008,0070 - Manufacturer tag)
    if dcm_data["Manufacturer"].value.lower().find('philips') > 0:
        scanner_type = 'philips'
    elif dcm_data["Manufacturer"].value.lower().find('siemens'):
        scanner_type = 'siemens'
    elif dcm_data["Manufacturer"].value.lower().find('ge'):
        scanner_type = 'ge'
    else:
        scanner_type = 'other'

    for tag_name in lut_tag:
        try:
            tag_value = str(dcm_data[tag_name].value)
            if isinstance(tag_value, list):
                tag_value = tag_value[0]
            if len(tag_value.split('\n')) > 2:
                logger.warning(f"Multiple values found for tag '{tag_name}' in file {str(dcm_path)}")
                for line in tag_value.split('\n'):
                    if (not line == '') and (line.find('(no value available)') == -1):
                        tag_value = line
                        break
                if tag_value.find('\n'):
                    tag_value = 'Unknown_' + tag_name
            # Initial replacing of bad characters into useful partitions
            tag_value = re.sub(r"[ ./']", '-', tag_value)
            tag_value = re.sub(r'\\', '_', tag_value)
            tag_value = re.sub(r'\*', 's', tag_value)
            tag_value = re.sub(r'\?', 'q', tag_value)
            tag_value = re.sub(r'\W', '', tag_value)
        except KeyError:
            logger.warning(f"WARNING: Missing DICOM tag '{tag_name}' in file {str(dcm_path)}")
            tag_value = 'Unknown_' + tag_name

        lut_value[tag_name] = tag_value

    dir_base = f"{lut_value['StudyInstanceUID']}_{lut_value['StudyDate']}_{lut_value['PatientName']}"

    if not lut_value["SeriesNumber"].startswith("Unknown"):
        lut_value["SeriesNumber"] = f'{int(lut_value["SeriesNumber"]):03d}'
    dir_series = f'{lut_value["SeriesNumber"]}-{lut_value["SeriesDescription"]}-UID_{lut_value["SeriesInstanceUID"]}'

    full_out = options.dir_output / dir_base / dir_series / f'{options.fname_dcm}.{options.extension}'
    full_out.parent.mkdir(parents=True, exist_ok=True)

    if full_out.exists() and not options.clobber:
        full_out = options.dir_output / dir_base / dir_series / f'{options.fname_dcm}A.{options.extension}'

    if options.move:
        dcm_path.rename(full_out)
    else:
        shutil.copyfile(dcm_path, full_out)
