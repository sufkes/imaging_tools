#!/usr/bin/env bash

module load dcmtk

scan_list="$1"
scan_type="$2"

## Convert the scans in each list to minc
base_dir="$PWD"

cat "$scan_list" | while read dd
do
    studyname="$(echo "$dd" | rev | cut -d / -f 2 | rev)"
    id_visit="$(echo "$studyname" | cut -d _ -f 2-3)"

    date="$(echo "$studyname" | cut -d _ -f 1)"
    seq_name="$(find "$dd" -type f | head -1 | while read ff; do dcmdump "$ff" | grep SeriesDescription | head -1 | cut -d [ -f 2 | cut -d ] -f 1; done)"
#    seq_name"$(basename "$dd" | cut -d _ f)"

    seq_name="$(basename "$dd" | while read sn; do python -c "print '""$sn""'.split('-UID')[0]"; done)"
    out_dir="$(echo "$id_visit" | cut -c 1-4)"_minc
    fname="$(echo "$id_visit" | tr '[:upper:]' '[:lower:]')"_"$date"_"$(echo "$seq_name" | tr " " _)"_"$scan_type"

    qsub <<EOF
#PBS -N $fname
#PBS -l mem=2g,vmem=2g
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
module load minc-toolkit
source /hpf/tools/centos6/minc-toolkit/1.0.05/minc-toolkit-config.sh
cd $base_dir
dcm2mnc "$dd" "$out_dir" -fname "$fname" -dname ''
EOF
done

 
