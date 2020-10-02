#!/usr/bin/env bash

usage() {
echo "Generate an FSL .mat and .con file for the specified columns numbers, and patients with the motion threshold code (025mm, 05mm, 1mm, 2mm, 4mm). This is designed for use in the MagNUM fMRI project, where subject motion was rated as <0.25mm, <0.5mm, <1mm etc.

Usage:
    makeMat.sh <data CSV> <colum numbers (e.g. 1-2)> <motion thresh code> <quantitative or categorical ('q' or 'c')> <output prefix>"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ ! -z "$6" ]
then
    usage
    exit 1
fi

module load fsl/6.0.0 &> /dev/null
source "$FSLDIR"/etc/fslconf/fsl.sh

data_csv="$1"
columns="$2"
motion_thresh="$3"
quant_or_cat="$4"
out_prefix="$5"
out_name="$out_prefix"_"$motion_thresh"

# Generate the .mat file.
tmpfile=$(mktemp)
cat "$data_csv" | grep "$motion_thresh" | cut -d "," -f "$columns" | tr "," '\t' > "$tmpfile"
Text2Vest "$tmpfile" "$out_name".mat
rm -f "$tmpfile"

# Generate the .con file.
tmpfile=$(mktemp)
if [ "$quant_or_cat" == "c" ]
then
    echo -e "1 -1\n-1 1" > "$tmpfile"
    Text2Vest "$tmpfile" "$out_name".con
elif [ "$quant_or_cat" == "q" ]
then
    echo -e "1\n-1" > "$tmpfile"
    Text2Vest "$tmpfile" "$out_name".con
fi
rm -f "$tmpfile"
