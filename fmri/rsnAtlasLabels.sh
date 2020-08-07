#!/bin/bash

# Written for FSL 6

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
    for label_num in $(seq 1 $num_labels)
    do
	label_mask=./${label_mask_prefix}${label_num}

	# Mask the IC using the current label.
	ic_masked=./ic${ic_num}_label${label_num}
	fslmaths "$ic_probmap" -mas "$label_mask" "$ic_masked"

	# Compute mean p of masked IC
	echo "$label_num" "$(fslstats "$ic_masked" -M)" "$(cat "$atlas_csv_path" | head -n "$label_num" | tail -n 1 | cut -d "," -f 2)" >> "$ic_report" # -M gets mean of non-zero voxels.
    done
done



