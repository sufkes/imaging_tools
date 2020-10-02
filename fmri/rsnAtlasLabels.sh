#!/bin/bash

# Written for FSL 6

usage() {
echo "Report the anatomical regions in which MELODIC ICAs appear.

Usage:

   rnsAtlasLabels.sh <MELODIC_DIR> <PATH_TO_ATLAS_NIFTI_FILE> <PATH_TO_ATLAS_LABEL_CSV>

MELODIC_DIR is the output directory of FSL MELODIC. It must contain a stats folder with probmap images.
PATH_TO_ATLAS_NIFTI is the path to a NIFTI file containing atlas labels. Labels must be integers. This image must be aligned with the images used in MELODIC
PATH_TO_ATLAS_LABEL_CSV is the path a CSV file which contains the label numbers in the first column, and the descriptions of the ROIs in the second column"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]
then
    usage
    exit 1
fi

melodic_dir="$1" # path to MELODIC output directory. Must contain a stats subdirectory with probmap images.
atlas_path="$2" # path to atlas label image (in same image space as probmap images) with integer labels > 1
atlas_csv_path="$3" # path to CSV with first column being the atlas label value, and the second column describing the ROI.

# Check that inputs exist
test -d "$melodic_dir" || { echo "MELODIC directory does not exist."; exit 1; }
test -f "$atlas_path" || { echo "Path to atlas does not exist."; exit 1; }
test -f "$atlas_csv_path" || { echo "Path to atlas CSV does not exist."; exit 1; }

# Guess the number of labels in the atlas label image.
num_labels="$(fslstats "$atlas_path" -R | awk '{print $2}' | cut -d . -f 1)" && echo "Assuming "$num_labels" labels exist in atlas label image."

# Generate a mask of each atlas label
label_mask_prefix=label_
for ii in $(seq 1 "$num_labels")
do
    lower=$(($ii - 1)).9
    upper=${ii}.1
    out_path=./${label_mask_prefix}${ii}
    fslmaths "$atlas_path" -thr $lower -uthr $upper "$out_path"
done

## Compute the average probabilty of "activation" within each atlas label for each MELODIC IC.
# Guess number of ICs
num_ics="$(ls "$melodic_dir"/stats/*probmap* | wc -l | awk '{print $1}')" && echo "Assuming "$num_ics" ICs."

for ic_num in $(seq 1 $num_ics)
do
    ic_report=./ic${ic_num}.txt
    
    ic_probmap="$melodic_dir"/stats/probmap_${ic_num}.nii.gz
    ic_threshz="$melodic_dir"/stats/thresh_zstat${ic_num}.nii.gz
    for label_num in $(seq 1 $num_labels)
    do
	label_mask=./${label_mask_prefix}${label_num}

	# Mask the IC using the current label.
	ic_probmap_masked=./ic_probmap_${ic_num}_label${label_num}
	fslmaths "$ic_probmap" -mas "$label_mask" "$ic_probmap_masked"

	ic_threshz_masked=./ic_threshz_${ic_num}_label${label_num}
	fslmaths "$ic_threshz" -mas "$label_mask" "$ic_threshz_masked"
	
	# Get the number of pixels of the IC in the current region for which p>0.5
	vol_p_ge_50="$(fslstats "$ic_probmap_masked" -l 0.50 -V | awk '{print $1}')"

	if [ $vol_p_ge_50 -ge 1 ] # only get more info if the current label region contains significant activiation.
	then
	
	    # Get the mean value of p within the label region.
	    mean_p="$(fslstats "$ic_probmap_masked" -M | awk '{print $1}')"

	    # Get the max value of p within the label region.
	    max_abs_z="$(fslstats "$ic_threshz_masked" -a -R | awk '{print $2}')"
	    
	    # Get the description of the current label region.
	    label_desc="$(cat "$atlas_csv_path" | head -n "$label_num" | tail -n 1 | cut -d "," -f 2)"
	
	    # Add the stats to the report.
	    echo "$label_num" "$max_abs_z" "$mean_p" "$label_desc" >> "$ic_report" # -M gets mean of non-zero voxels.
	fi

	# Remove the masked IC.
	rm "$ic_probmap_masked".nii.gz
	rm "$ic_threshz_masked".nii.gz
    done
# Sort the IC report file according to the mean p values.
sort -o "$ic_report" -k 2 -r "$ic_report"
done

# Remove the label mask files.
for ii in $(seq 1 "$num_labels")
do
    out_path=./${label_mask_prefix}${ii}
    rm "$out_path".nii.gz
done

echo "Run combineAtlasLabels.sh to combine results in a single CSV file."
