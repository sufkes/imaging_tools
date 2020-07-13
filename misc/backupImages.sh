#!/bin/bash

# Parse arguments.
usage() {
    help="Usage:\n\tbackupImages.sh COHORT_DIR\n\nCOHORT_DIR is a directory containing DICOM Study subdirectories and possibly some files. Each subdirectory in COHORT_DIR will be compressed and uploaded to the iRods Archive. Files in COHORT_DIR will be uploaded without compression.\n\nNote: This script must be run in an interactive session in which the 'irods_client' module is available. Currently, the only way to do this is to log into hpf.ccm.sickkids.ca, and start an interactive session on a centos7 node using 'qsub -q centos7 -I'."
    echo -e $help
    exit 1
    }

if [ -z "$1" ]
then
    usage
else
    cohort_dir="$1"
fi

# Load iRods module.
module_error() {
    msg="Error: Failed to load 'irods_client' module.\n\nNote: This script must be run in an interactive session in which the 'irods_client' module is available. Currently, the only way to do this is to log into hpf.ccm.sickkids.ca, and start an interactive session on a centos7 node using 'qsub -q centos7 -I'."
    echo -e $msg
    exit 2
    }
module load irods_client 2>&1 | grep -q ERROR && module_error # module load returns 0 if module does not exist, so grep for the error message instead of looking at the exit code.

# Check that cohort directory exists and is not a sumbolic link.
if [ -d "$cohort_dir" ] && [ ! -L "$cohort_dir" ]
then
    :
    #echo "Backing up: ${cohort_dir}"
else
    echo "Error: Input cohort directory is not a directory: ${cohort_dir}"
    exit 3
fi

# Determine destination on Archive using absolute path of cohort directory. (Replace "/hpf/largeprojects/smiller/images" with /"resarchivezone/smiller_grp_images_backup")
cohort_dir_abs_path="$(readlink -f "$cohort_dir")"
archive_dir="$(echo "$cohort_dir_abs_path" | sed 's/\/hpf\/largeprojects\/smiller\/images/\/resarchivezone\/smiller_grp\/images_backup/')"

# Check whether destination on archive already exists. This should check if it is a collection (instead of a dataObj), but I don't know how to do that nicely, and it doesn't seem likely to cause problems.
confirm() {
    read -p "Enter (Y/y) to continue: " -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]
    then
	exit 4
    fi
    }
ils "$archive_dir" &> /dev/null
if [ $? -eq 0 ]
then
    echo "Destination already exists on Archive. This data will be overwritten: ${archive_dir}"
else
    echo "Destination does not exist on Archive. It will be created: ${archive_dir}"
fi
confirm # wait for user confirmation.

# Make temporary collection in Archive to save files to as they are compressed and uploaded. Check if it already exists. Warn about creation.
ils /resarchivezone/smiller_grp/temp &> /dev/null || imkdir /resarchivezone/smiller_grp/temp
temp_archive=/resarchivezone/smiller_grp/temp/"$(basename "$cohort_dir_abs_path")"
ils "$temp_archive" &> /dev/null
if [ $? -eq 0 ]
then
    echo "Error: Temporary collection already exists in Archive: ${temp_archive}"
    echo "Remove or rename this collection and restart."
    exit 5
else
    echo "Creating temporary collection in Archive: ${temp_archive}"
    echo "If backup fails, you should delete this temporary collection manually."
fi

# Add HPF data to temporary Archive collection.
imkdir "$temp_archive"

compressAndUpload() {
    in="$1" # input file or directory
    out="$2" # output collection in Archive

    in_basename="$(basename "$in")"
    in_dirname="$(dirname "$in")"
    
    if [ ! -L "$in" ] # Don't do anything with symbolic links.
    then
	qsub <<EOF
#!/bin/bash
#PBS -N compressAndUpload
#PBS -l walltime=12:00:00
#PBS -l mem=8g,vmem=8g
#PBS -l nodes=1:ppn=1
#PBS -j oe

module load irods_client || exit 1

# Determine whether input is a file or a directory. Assume that file is not a symbolic link.
if [ -f "$in" ]
then 
    # Add file to temporary collection in Archive immediately.
    iput "$in" "$out" || exit 1
elif [ -d "$in" ]
then
    # Create temporary file to store the compressed directory.
    temp_file="\$(mktemp)" || exit 1

    # Compress the input directory.
    tar -zcf "\$temp_file" -C "$in_dirname" "$in_basename" || exit 1

    # Add compressed directory to temporary collection in Archive.
    iput "\$temp_file" "$out"/"$(basename "$in")".tar.gz || exit 1

    # Delete the temporary file.
    rm "\$temp_file"
fi
exit 0
EOF
    fi
    
}

dependency_jobs=""
while read input
do
    job_id="$(compressAndUpload "$input" "$temp_archive")"

    # Add JOB ID to list of dependency jobs.
    if [ ! -z "$job_id" ]
    then
	dependency_jobs="$dependency_jobs":"$job_id"
    fi
done < <(find "$cohort_dir_abs_path" -mindepth 1 -maxdepth 1 | sort)

# Once all files have been added to the temporary archive location, remove the previous backup, and replace it with the new one.
qsub -W depend=afterok"$dependency_jobs" <<EOF > /dev/null
#!/bin/bash
#PBS -N "$(basename "$cohort_dir_abs_path")"
#PBS -l walltime=24:00:00
#PBS -l mem=2g,vmem=2g
#PBS -l nodes=1:ppn=1
#PBS -m a
#PBS -j oe

module load irods_client || exit 1

# Delete the old backup.
ils "$archive_dir" &> /dev/null
if [ $? -eq 0 ]
then
    irm -r "$archive_dir" || exit 1
fi

# Move the new backup to the proper location.
imv "$temp_archive" "$archive_dir" || exit 1

EOF


