#!/usr/bin/env bash

module load dcmtk

t1=./t1_scans.txt
t2=./t2_scans.txt

## Email from Jessie 2020-05-15:
# The sequences for T1-weighted images are normally contain the following key words:
# - T1-Ax-3D-FLASH
# - mpr or MPR (and the name doesn't contain "T2")
# - Sag-fl3D1r

# The sequences for T2-weighted images are normally contain the following key words:
# - T2-Ax-2D-TSE
# - T2 & mpr or T2 & MPR

## T1 scans
# T1-Ax-3D-FLASH -> This will include 2 scans missing 3D in their name: 20161219_SK040019_V01/101-T1-flash-cor-mpr-1mm and 20161219_SK040019_V01/100-T1-flash-ax-mpr)
find -L . -mindepth 3 -maxdepth 3 | grep -i t1 | grep -i flash > "$t1"

# mpr or MPR (and the name doesn't contain "T2")
find -L . -mindepth 3 -maxdepth 3 | grep -i mpr | grep -vi t2 >> "$t1"

# Sag-fl3D1r
find -L . -mindepth 3 -maxdepth 3 | grep -i fl3D1r >> "$t1"


## T2 scans
# T2-Ax-2D-TSE
find -L . -mindepth 3 -maxdepth 3 | grep -i t2 | grep -i tse > "$t2" # (T2 or t2) and (TSE or tse)

# T2 & mpr or T2 & MPR
find -L . -mindepth 3 -maxdepth 3 | grep -i t2 | grep -i mpr >> "$t2"


## Clean up lists.
t1_temp="$(mktemp)"
cat "$t1" | sort | uniq > "$t1_temp"
cat "$t1_temp" > "$t1"
rm "$t1_temp"
t2_temp="$(mktemp)"
cat "$t2" | sort | uniq > "$t2_temp"
cat "$t2_temp" > "$t2"
rm "$t2_temp"

## Look for scans that lie in both lists.
cat "$t1" | while read dd;
do
    cat "$t2" | grep -q "dd"
    if [ $? = 0 ]
    then
	echo 'Scan in both lists:' "$dd"
    fi
done
    
    
