#!/bin/bash

#PBS -N mr_diamond
#PBS -l nodes=1:ppn=16
#PBS -l mem=32g,vmem=32g
#PBS walltime=02:00:00
#PBS -joe
#PBS -m e

# Usage
usage () {
	echo "Usage:"
	echo '    qsub mrdiamond.sh -F "abs_dir_path subject_id visit_id"'
}
# To pass arguments into script submitted using qsub:
# qsub mr_sort.sh -F "in_dir subject_id visit_id"

# Check that arguments are supplied
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]
then
	echo "Missing arguments"
	usage # print usage information.
	exit 1
fi

module load dcmtk
module load parallel
export PATH=$PATH:/hpf/largeprojects/smiller/users/sufkes/sufkes_imaging/dicom
export PATH=$PATH:/hpf/largeprojects/smiller/users/sufkes/sufkes_imaging/misc

NP=$(wc -l $PBS_NODEFILE | awk '{print $1}') # Should return the number of physical CPU cores.

# Change to directory containing files to sort.
in_dir="$1"
cd "$in_dir"

# Get subject ID and Visit ID
subject="$2"
visit="$3"

# Count files
num_files_start="$(find . -type f | wc -l)"
echo "Attempting to sort ""$num_files_start"" files."

# Make directories
echo "Making directories for sorted files."
mkdir raw
mkdir dicom
mkdir rda

# Move all of the unsorted data to the raw folder.
echo "Moving unsorted files."
find . -type f | while read ff; do mv -n "$ff" raw; done

# Sort DICOMs
echo "Sorting DICOM files."
find raw -type f | parallel -j "$NP" processDicom.py "{}" dicom -n "$subject" -i "$visit"

# Anonymize DICOMs
echo "Anonymizing DICOM files."
find dicom -type f | parallel -j "$NP" anonDicom.py "{}" -n "$subject" -m

# Copy RDA
echo "Copying RDA files."
find raw -type f | grep ".rda" | while read ff; do cp -n "$ff" rda; done

# Anonymize RDA
echo "Anonymizing RDA files."
find rda -type f | while read ff; do anonRDA.py "$ff" -n "$subject"; done

# Rename RDA folder
echo "Renaming RDA folder and files."
new_name="$(find dicom -maxdepth 1 -mindepth 1 -type d | sort | head -1 | xargs basename | cut -d _ -f 1-3)"
mkdir rda/"$new_name"
mv -n rda/*rda rda/"$new_name"

# Change ownership & permissions
echo "Changing ownership and permissions."
chown -R sufkes:smiller dicom rda
chmod -R 550 dicom/* rda/*

# Check final number
num_files_end="$(find dicom rda -type f | wc -l)"
#echo "Number of sorted files:" "$num_files_end"
if [ ! "$num_files_start" = "$num_files_end" ]
then
	echo "Warning: Number of sorted files differs from number of original files."
else
	echo "Finished sorting ""$num_files_end"" files."
fi
