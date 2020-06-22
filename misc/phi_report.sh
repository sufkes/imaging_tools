#!/bin/bash

#PBS -l walltime=12:00:00
#PBS -l mem=8g,vmem=8g
#PBS -joe /home/sufkes/logs/phi_reports
#PBS -me 

check_dir="$1"
out_path="$check_dir"/../"$(basename "$check_dir")"_phi_report_deep.txt
cd "$check_dir"
ls -d */ | tr -d / | while read dd; do echo "$dd"; echo "DICOM tags:"; find "$dd" -type f | while read ff; do anonDicom.py -p "$ff"; done | sort | uniq; echo; echo "Possible requisition forms:"; find "$dd" -type d | grep -Ei 'req|paper|999|1000'; echo "--------------------------------------------------------------------"; done > "$out_path"
