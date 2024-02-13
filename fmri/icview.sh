#!/bin/bash

usage() {
echo """View independent components from MELODIC for a single subject.

Usage:
    d <study name> [ICs to view]

Example 1: View all ICs for a single subject:
./d 20210420_BC0266_V08

Example 2: View the fifth component for a single subject:
./d 20210420_BC0266_V08 5

Example 3: View components 1, 4, and 5 for a single subject:
./d 20210420_BC0266_V08 1 4 5"""
}


# Parse arguments
study="$1"
background_path="${study}/${study}_ref_vol.nii.gz"

# FSLEyes display settings
pos_color='red-yellow'
neg_color='blue-lightblue'

if [ $# -eq 0 ]; then
    usage
    exit
elif [ $# -eq 1 ]; then
    # Get paths to all positive components.
    num_ics=$(find "${study}/stats" -name 'pos_*nii.gz' | wc -l)
    pos_paths=() # array storing the paths to the positive parts of the ICs, in reverse order (to make them in the correct order in fsleyes).
    for index in $(seq "$num_ics" 1); do
	pos_path="${study}/stats/pos_thresh_zstat${index}.nii.gz"
	pos_paths+=("$pos_path")
    done
    hide=true # FSLeyes argument to hide image by default.
elif [ $# -eq 2 ]; then
    # Get path to the requested IC.
    ic="$2"
    pos_paths=("${study}/stats/pos_thresh_zstat${ic}.nii.gz")
    hide=false
else
    # The paths to the requested ICs.
    args=("$@")
    pos_paths=()
    for arg_num in $(seq $(($# - 1)) 1); do # loop over the arguments, starting with the second one, in reverse order.
	ic="${args[$arg_num]}"
	pos_paths+=("${study}/stats/pos_thresh_zstat${ic}.nii.gz")
    done
    hide=true
fi

# Make the FSLeyes command arg-by-arg.
cmd=(fsleyes "$background_path")
for pos_path in "${pos_paths[@]}"; do
    neg_path="$(dirname "$pos_path")/$(basename "$pos_path" | sed "s/^pos/neg/g")"

    pos_min=$(fslstats "$pos_path" -P 0) # 0th percentile of nonzero voxels (i.e. min of nonzero voxels).
    pos_max=$(fslstats "$pos_path" -P 100)
    neg_min=$(fslstats "$neg_path" -P 0)
    neg_max=$(fslstats "$neg_path" -P 100)

    cmd+=("$neg_path" -cm "$neg_color" -dr "$neg_min" "$neg_max")
    if $hide; then cmd+=("-d"); fi    
    cmd+=("$pos_path" -cm "$pos_color" -dr "$pos_min" "$pos_max")
    if $hide; then cmd+=("-d"); fi
done

echo "${cmd[@]}"
"${cmd[@]}"
