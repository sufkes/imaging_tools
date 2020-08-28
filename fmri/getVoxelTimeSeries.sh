#!/usr/bin/env bash

usage() {
    echo "Extract fMRI time series for every voxel in an ROI"
    echo
    echo "Usage:"
    echo "    getVoxelTimeSeries.sh <fMRI_file_list> <ROI_mask_list> <out_dir>"
    exit 1
}

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]
then
    usage
fi


bold_list="$1" # text file containing the paths to the BOLD images
mask_list="$2" # text file containing the paths to the ROI masks
out_dir="$3" # path to directory where time series will be saved.

module load fsl/6.0.0
source "$FSLDIR"/etc/fslconf/fsl.sh

while read bold
do
    while read mask
    do
	out_path="$out_dir"/"$(basename "$bold" | cut -d . -f 1)"_"$(basename "$mask" | cut -d _ -f 2)".txt
	fslmeants -i "$bold" -m "$mask" -o "$out_path" --showall 
    done < <(cat "$mask_list")
done < <(cat "$bold_list")
