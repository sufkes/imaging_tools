#!/usr/bin/env python

# Standard modules
import os, sys
from collections import OrderedDict
import argparse
import subprocess

# Non-standard modules
import pandas as pd

def run_cmd(sys_cmd, verbose=False):
    """Run a terminal command"""
    if verbose:
        print sys_cmd
    p = subprocess.Popen(sys_cmd, stdout = subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = p.communicate()
    return output, errors

def isDicom(in_path):
    """Check if file is DICOM using the 'file' command."""
    cmd_file = 'file "{in_path}"'.format(in_path=in_path)
    output, errors = run_cmd(cmd_file, 0)
    if ("DICOM medical imaging data" in output):
        return True
    else:
        return False

def buildPathList(in_dirs): # get paths of all files in in_dirs
    """Build list of absolute paths of all DICOM files in specified directories."""
    pathList = []
    for in_dir in in_dirs:
        for dirpath, dirnames, filenames in os.walk(in_dir):
            for filename in filenames:
                rel_path = os.path.join(dirpath, filename)

                # Check if file is DICOM.
                if (not isDicom(rel_path)):
                    print "Warning: File appears to be non-DICOM or corrupted. Excluding from report: '{file}'".format(file=rel_path)
                    continue
                
                # Convert relative path to absolute path.
                abs_path = os.path.abspath(rel_path)

                # Add path to list.
                pathList.append(abs_path)

    # Sort path list
    pathList = sorted(pathList, key=lambda p: (os.path.sep not in p, p))
    return pathList

def getDicomTag(in_path, tag_name):
    """Get the value stored in a DICOM tag using 'dcmdump'.
    """

    # Buid dict of tag name -> number mappings.
    tag_dict = {}
    tag_dict["StudyDate"] = "0008,0020"
    tag_dict["StudyTime"] = "0008,0030"
    tag_dict["SeriesNum"] = "0020,0011"
    tag_dict["StudyDescription"] = "0008,1030"
    tag_dict["ProtocolName"] = "0018,1030"
    tag_dict["InstanceNumber"] = "0020,0013"
    tag_dict["SeriesDescription"] = "0008,103e"
    tag_dict["PatientName"] = "0010,0010"
    tag_dict["StudyInstanceUID"] = "0020,000d"
    tag_dict["SeriesInstanceUID"] = "0020,000e"
    tag_dict["MRAcquisitionType"] = "0018,0023"
    tag_dict["Manufacturer"] = "0008,0070" # seems always present, not always consistently named
    tag_dict["ManufacturerModelName"] = "0008,1090" # seems always present, not always consistently named
    tag_dict["SequenceName"] = "0018,0024" # not always present
    tag_dict["ImageType"] = "0008,0008" # e.g. raw DTI: [ORIGINAL\PRIMARY\DIFFUSION\NONE\ND\MOSAIC], derived DTI: [DERIVED\PRIMARY\DIFFUSION\ADC\ND] 
    
    # Convert tag name to a tag number:
    tag_number = tag_dict[tag_name]

    # Get the tag value using dcmdump
    cmd_dcmdump = 'dcmdump "{in_path}" +P "{tag_number}" -s +L | cut -d [ -f 2- | rev | cut -d ] -f 2- | rev'.format(in_path=in_path, tag_number=tag_number) # Get complete value in first instance of tag. Grab everything lying between the first "[" and the last "]".
    out, err = run_cmd(cmd_dcmdump)
    out = out.rstrip("\n") # the output of the command has a newline appended to it.

    # Check if tag has (no value available), or tag does not appear.
    if (out == "") or ("(no value available)" in out):
        value = None
    else:
        value = out
    return value

series_types = ["3DT1",
                "2DT2",
                "DTI",
                "DTI_derived",
                "DWI",
                "fMRI",
                "unknown"]

class Series(object): # information about the DICOM Series.
    def __init__(self, SeriesInstanceUID):
        self.SeriesInstanceUID = SeriesInstanceUID # Can only have a single value by construction.

        self.series_dir = set() # hopefully will only have one value.
        self.paths = [] # paths to files in this series.

        self.SeriesNum = set() # hopefully one value
        self.SeriesDescription = set() # hopefully one value
        self.ProtocolName = set() # hopefully one value
        self.SequenceName = set() # hopefully one value
        self.ImageType = set() # hopefully one value        
        self.MRAcquisitionType = set() # hopefully one value
        self.Manufacturer = set() # hopefully one value
        self.ManufacturerModelName = set() # hopefully one value
        self.InstanceNumbers = [] # likely many values

        self.num_files = 0
        self.maxInstanceNumber = None

        ## Need to handle this better. These need to be publicly available, and more easily modifiable.
        self.type_dict = OrderedDict([("3DT1", self.__is_3DT1),
                                      ("2DT2", self.__is_2DT2),
                                      ("DTI", self.__is_DTI),
                                      ("DTI_derived", self.__is_DTI_derived),
                                      ("DWI", self.__is_DWI),
                                      ("fMRI", self.__is_fMRI),
                                      ("unknown", self.__is_unknown)])
        assert set(self.type_dict.keys()) == set(series_types)
        self.series_type = set() # hopefully one value

        self.all_files_examined = False # Whether all files in the series have been added to this object.
        self.num_files_equals_max_instance = False # Whether the number of files is equal to the maximum instance number.
        
    def addDir(self, file_path):
        # Assume that Series directory is the directory containing the current file.
        series_dir = os.path.dirname(file_path)
        self.series_dir.add(series_dir)
        
    def addPath(self, file_path):
        self.paths.append(file_path)

    def addSeriesNum(self, file_path):
        SeriesNum = int(getDicomTag(file_path, "SeriesNum"))
        self.SeriesNum.add(SeriesNum)
        
    def addSeriesDescription(self, file_path):
        SeriesDescription = getDicomTag(file_path, "SeriesDescription")
        if (not SeriesDescription is None):
            self.SeriesDescription.add(SeriesDescription)

    def addProtocolName(self, file_path):
        ProtocolName = getDicomTag(file_path, "ProtocolName")
        if (not ProtocolName is None):
            self.ProtocolName.add(ProtocolName)

    def addSequenceName(self, file_path):
        SequenceName = getDicomTag(file_path, "SequenceName")
        if (not SequenceName is None):
            self.SequenceName.add(SequenceName)

    def addImageType(self, file_path):
        ImageType = getDicomTag(file_path, "ImageType")
        if (not ImageType is None):
            self.ImageType.add(ImageType)

    def addMRAcquisitionType(self, file_path):
        MRAcquisitionType = getDicomTag(file_path, "MRAcquisitionType")
        if (not MRAcquisitionType is None):
            self.MRAcquisitionType.add(MRAcquisitionType)

    def addManufacturer(self, file_path):
        Manufacturer = getDicomTag(file_path, "Manufacturer")
        if (not Manufacturer is None):
            self.Manufacturer.add(Manufacturer)

    def addManufacturerModelName(self, file_path):
        ManufacturerModelName = getDicomTag(file_path, "ManufacturerModelName")
        if (not ManufacturerModelName is None):
            self.ManufacturerModelName.add(ManufacturerModelName)

    def addInstanceNumber(self, file_path):
        InstanceNumber = int(getDicomTag(file_path, "InstanceNumber"))
        if (not InstanceNumber is None):
            self.InstanceNumbers.append(InstanceNumber)
    
    def addFile(self, file_path):
        self.addDir(file_path)
        self.addPath(file_path)
        self.addSeriesNum(file_path)
        self.addSeriesDescription(file_path)
        self.addProtocolName(file_path)
        self.addSequenceName(file_path)
        self.addImageType(file_path)
        self.addMRAcquisitionType(file_path)
        self.addManufacturer(file_path)
        self.addManufacturerModelName(file_path)
        self.addInstanceNumber(file_path)

    def setAllFilesExamined(self, value=True):
        self.all_files_examined = value
        
    def checkAllFilesExamined(self):
        if (not self.all_files_examined):
            raise Exception("Cannot do this if Series.all_files_examined is False")
        
    def setSummaryValues(self):
        """Calucate measures for this series which require that data has been added for all files in this series."""
        self.checkAllFilesExamined() # ensure that all files have been examined.
        self.num_files = len(self.paths)
        self.maxInstanceNumber = max(self.InstanceNumbers)

    def setQualityMetrics(self):
        self.checkAllFilesExamined() # ensure that all files have been examined.
        self.num_files_equals_max_instance = (self.num_files == self.maxInstanceNumber) # weak check to see whether all files in a Series are present.
        
    def __is_3DT1(self):
        # Jessie says: The sequences for T1-weighted images are normally contain the following key words:
        # - T1-Ax-3D-FLASH
        # - mpr or MPR (and the name doesn't contain "T2", but could contain "rpt2" or "repeat2")
        # - Sag-fl3D1r

        ## Get the SeriesDescription and MRAcquisitionType.
        # Convert uppercase to lowercase.
        # If more than one unique values were found, concatenate them.
        # Remove strings "rpt" or "repeat" as they add no useful information, and can cause problems (e.g. if SeriesDescription says "RPT2" to indicate a second repetition of a sequence)
        desc = " ".join(self.SeriesDescription).lower().replace("rpt", "").replace("repeat", "")
        mracq = " ".join(self.MRAcquisitionType).lower()

        # Does it look like T1?
        if ("t1" in desc) or (("mpr" in desc) and (not "t2" in desc)) or ("fl3d1r" in desc):
            # Does it look like 3D?
            if ("3d" in mracq):
                return True
        return False

    def __is_2DT2(self):
        # Jessie says: The sequences for T2-weighted images are normally contain the following key words:
        # T2-Ax-2D-TSE
        # T2 & mpr
        
        ## Get the SeriesDescription and MRAcquisitionType.
        desc = " ".join(self.SeriesDescription).lower().replace("rpt", "").replace("repeat", "")
        mracq = " ".join(self.MRAcquisitionType).lower()
        # Does it look like T2?
        if ("t2" in desc) and (("tse" in desc) or ("mpr" in desc)):
            # Does it look like 2D?
            if ("2d" in mracq):
                return True
        return False

    def __is_GRE(self):
        ## Get the SeriesDescription and ImageType
        desc = " ".join(self.SeriesDescription).lower().replace("rpt", "").replace("repeat", "")
        if (("gre" in desc) or ("field" in desc) or ("map" in desc)):
            return True
    
    def __is_DTI(self):
        ## Get the SeriesDescription and ImageType
        desc = " ".join(self.SeriesDescription).lower().replace("rpt", "").replace("repeat", "")
        imtype = " ".join(self.ImageType).lower()
        # Does it look like DTI?
        if ("dti" in desc):
            # Does it look like a GRE field mapping?
            if (not self.__is_GRE()):
                # Does it look like raw DTI, as opposed to a DTI-derived image (e.g. ADC, FA)?
                if ("original" in imtype):
                    return True
        return False

    def __is_DTI_derived(self):
        ## Get the SeriesDescription and ImageType
        desc = " ".join(self.SeriesDescription).lower().replace("rpt", "").replace("repeat", "")
        imtype = " ".join(self.ImageType).lower()
        # Does it look like DTI
        if ("dti" in desc):
            # Does it look like a GRE field mapping?
            if (not self.__is_GRE()):
                # Does it look like a DTI-derived image, as opposed to raw DTI?
                if ("derived" in imtype):
                    return True
        return False

    def __is_DWI(self):
        return False

    def __is_fMRI(self):
        ## Get the SeriesDescription and MRAcquisitionType.
        desc = " ".join(self.SeriesDescription).lower().replace("rpt", "").replace("repeat", "")
        mracq = " ".join(self.MRAcquisitionType).lower()
        # Does it look like fMRI, and not a GRE field map?
        if (("fmri" in desc) or ("fcmri" in desc) or ("resting" in desc) or ("bold" in desc)):
            # Does it look like a GRE field mapping?
            if (not self.__is_GRE()):
                # Does it look like 2D?
                if ("2d" in mracq):
                    return True
        return False

    def __is_unknown(self):
        # made these protected because __is_unknown should only be called in self.classify, after all the other classifications have been executed.
        if (self.series_type == set()):
            return True
    
    def classify(self):
        self.checkAllFilesExamined() # ensure that all files have been examined.
        for type_name, type_func in self.type_dict.iteritems():
            if type_func():
                self.series_type.add(type_name)

    def summarize(self):
        """Get information which requires that data has been extracted from all files in the series."""
        self.checkAllFilesExamined() # ensure that all files have been examined.

        # Set summary values (e.g. number of files)
        self.setSummaryValues()

        # Attempt to check quality of series
        self.setQualityMetrics()
        
        # Classify the series (as DTI, 3D T1 etc.).
        self.classify()

class Study(object): # information about DICOM Study (i.e. about the "scan")
    def __init__(self, StudyInstanceUID):
        self.StudyInstanceUID = StudyInstanceUID

        self.all_files_examined = False # whether all files in this Study have been examined.

        self.series = OrderedDict() # stores instances of Series object.
        self.study_dir = set() # stores all of the study directories identified.

        self.StudyDate = set() # hopefully one value
                
    def addDir(self, file_path):
        # Assume that the Study directory is one level above the directory containing the current file.
        study_dir = os.path.abspath(os.path.join(os.path.dirname(file_path), ".."))
        self.study_dir.add(study_dir)

    def addStudyDate(self, file_path):
        StudyDate = int(getDicomTag(file_path, "StudyDate"))
        self.StudyDate.add(StudyDate)
        
    def addFile(self, file_path):
        SeriesInstanceUID = getDicomTag(file_path, "SeriesInstanceUID")
        # Add information about the current file to this Study object (will often be redundant).
        self.addDir(file_path)
        self.addStudyDate(file_path)
        
        # Add information about current file to an instance of the Series object.
        if (not SeriesInstanceUID in self.series.keys()):
            # Create an instance of the Series object and add it to this study's dict of series.
            self.series[SeriesInstanceUID] = Series(SeriesInstanceUID)

        self.series[SeriesInstanceUID].addFile(file_path)

    def setAllFilesExamined(self, value=True):
        for SeriesInstanceUID, series in self.series.iteritems():
            series.setAllFilesExamined(value=value)
        self.all_files_examined = value

    def checkAllFilesExamined(self):
        if (not self.all_files_examined):
            raise Exception("Cannot do this if Series.all_files_examined is False")

    def summarizeSeriesInStudy(self):
        """Summarize information about the series in this study, once all the data has been extracted."""
        self.checkAllFilesExamined() # ensure that all files in this Study have been examined.
        for SeriesInstanceUID, series in self.series.iteritems():
            series.summarize()
            
def reportClassifications(studies):
    # Define helper for printing tabbed lines.
    def space(n, spacer="  "):
        return n*spacer
    
    for StudyInstanceUID, study in studies.iteritems():
        study_dir = list(study.study_dir)[0]
        print os.path.basename(study_dir)
        for series_type in series_types:
            print space(1)+series_type
            for SeriesInstanceUID, series in study.series.iteritems():
                if (series_type in series.series_type):
                    series_dir = list(series.series_dir)[0]
                    print space(2)+os.path.basename(series_dir)
            
def classifyDicom(in_dirs):
    # Get path to all DICOM files in in_dirs
    path_list = buildPathList(in_dirs)

    studies = OrderedDict()
    
    # Add information from each file.
    for path in path_list:
        # Get StudyInstanceUID
        StudyInstanceUID = getDicomTag(path, "StudyInstanceUID")

        # Add DICOM Study to Study dict if it is not yet there.
        if (not StudyInstanceUID in studies):
            # Initialize an instance of the Study object.
            studies[StudyInstanceUID] = Study(StudyInstanceUID)
            
        # Add information about this file to the Study.
        studies[StudyInstanceUID].addFile(path)

    for StudyInstanceUID, study in studies.iteritems():
        # Set the variables which record whether all files in the study/series have been examined. This is to ensure that measures which require knowledge of all files in the study/series to be known are not calculated prematurely.
        study.setAllFilesExamined() 

        # Calculate summary measures and classify the series using the information taken from the files.
        study.summarizeSeriesInStudy()

    # Report the classification results.
    reportClassifications(studies)
        
    # Print debug information.
    debug = False
    if debug:
        for StudyInstanceUID, study in studies.iteritems():
            print 80*"="
            # Print information about the Study.
            print "StudyInstanceUID:", study.StudyInstanceUID
            print "study_dir:", study.study_dir
            print "StudyDate:", study.StudyDate
            print "all_files_examined:", study.all_files_examined
            
            # Print information about the first few Series in the Study.
            for SeriesInstanceUID, series in study.series.iteritems():
                print 60*"-"
                print "series_dir:", series.series_dir
                #print "paths:", series.paths
                print "SeriesDescription:", series.SeriesDescription
                #print "ProtocolName:", series.ProtocolName
                #print "SequenceName:", series.SequenceName
                print "ImageType:", series.ImageType
                print "MRAcquisitionType:", series.MRAcquisitionType
                #print "Manufacturer:", series.Manufacturer
                #print "ManufacturerModelName:", series.ManufacturerModelName
                #print "InstanceNumbers:", series.InstanceNumbers
                print "series_type:", series.series_type
                #print "all_files_examined:", series.all_files_examined
                #print "num_files:", series.num_files
                #print "maxInstanceNumber:", series.maxInstanceNumber
                #print "num_files_equals_max_instance:", series.num_files_equals_max_instance
        print 

if (__name__ == '__main__'):
    # Create argument parser
    description = """Create a report on the types of Series present in DICOM Studies for a set of DICOM files. """
    parser = argparse.ArgumentParser(description=description)
    
    # Define positional arguments.
    parser.add_argument("in_dirs", help="path to directories containing DICOM files to be classified", type=str, nargs="+")
    
    # Define optional arguments.
#    parser.add_argument("-n", "--name", help="Subject ID to set PatientName tags to", type=str)
#    parser.add_argument("-b", "--backup", help="add an option for backup. Do not back up by default.")
#    parser.add_argument("-l", "--level", type=int, default=2, choices=[1,2,3], help="Set degree of anonymization (1: directly identifying information such as patient name, birth date, address, phone number. 2 (default): indirectly identifying information such as weight, age, physicians etc. 3: information about institution which performed scan, such as address, department etc.)")
#    parser.add_argument('-m', "--modify_pid", help="Change PatientID to specified Subject ID. Default: False", action="store_true")
#    parser.add_argument('-p', '--print-only', help='Print PHI-containing tags. Do not anonymize.', action='store_true')
#    parser.add_argument('-r', '--recursive', action='store_true', help='if in_path is a directory, find and anononymize all files in that directory.')
#    parser.add_argument("-r", "--records", help="list of records to export. Default: Export all records.", nargs="+", metavar=("ID_1", "ID_2"))

    # Parse arguments.
    args = parser.parse_args()

    # Classify the input DICOM files.
    classifyDicom(args.in_dirs)
