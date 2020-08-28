#!/usr/bin/env bash

usage() {
    echo "Split a single NIFTI file with n>0 integer labels into n images with labels 1 (ROI) and 0. Requires FSL."
    echo
    echo "Usage:"
    echo "    splitAtlasLabels.sh <NIFTI_file> [out_dir]"
}

test -z $1 && { usage; exit 1; }

in_path="$1"
if [ ! -z "$2" ]
then
    out_dir="$2"
else
    out_dir="$PWD"
fi

file_prefix="$(basename "$in_path" | cut -d . -f 1)"_label

echo "Reading label image: $in_path"
echo "Saving split label images to: $out_dir"

# Get maximum label value:
max="$(fslstats "$in_path" -R | awk '{print $2}' | cut -d . -f 1)"
if [ $max -lt 10 ] # zero pad the output file names.
then
    pad=1
elif [ $max -lt 100 ]
then
    pad=2
elif [ $max -lt 1000 ]
then
    pad=3
else
    pad=4
fi

for label in $(seq 1 $max)
do
    # grab voxels in range [label-0.1, label+0.1].
    thr="$(python -c "print(float("$label")-0.1)")" # lower threshold
    uthr="$(python -c "print(float("$label")+0.1)")" # upper threshold

    out_path="$out_dir"/"$file_prefix""$(printf "%0"$pad"d\n" $label)"
    
    fslmaths "$in_path" -thr "$thr" -uthr "$uthr" "$out_path"
done
