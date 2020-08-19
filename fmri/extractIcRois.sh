#!/usr/bin/env bash

melodic_dir="$1"
thresh="$2"

module load fsl/6.0.0
source "$FSLDIR"/etc/fslconf/fsl.sh

out_dir="$melodic_dir"/ica_rois
mkdir -p "$out_dir"

# For each ICA, get the probmap thresholded at some value, then binarize it to create an ROI.
find "$melodic_dir"/stats/ -name *probmap* | sort | while read ff
do
    nn="$out_dir"/"$(basename "$ff" | cut -d . -f 1)"_mask_thr"$(echo "$thresh" | tr -d '.')"
    fslmaths "$ff" -thr "$thresh" -bin "$nn"
done
