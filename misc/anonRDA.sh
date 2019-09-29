#!/usr/bin/env bash

# Anonymize RDA file header in place.
# Parameters:
#    ff : path to RDA file
#    nn (optional) : if specified, overwrite PatientName and PatientID with nn. If not specified, the PatientID is set to the value stored in the PatientName field of the header, and the PatientName is left unchanged.
# doesn't quite work; there is a bug in the actual modification part. Line with sed works but adds newline at end of file; line with perl doesn't work.

echo "This script doesn't quite work. Use anonRDA.py instead."
exit

ff="$1" # Path to RDA file
nn="$2" # Text to replace PatientName and PatientID with

change_tag() {
    local tag="$(grep -a "$1" "$ff" | cut -d " " -f 2- | cut -d $'\r' -f 1)"
    local tag_line="$(grep -a -n -m 1 "$1" "$ff" | cut -d ':' -f 1)"
    local str_to_search="$(grep -a -m 1 "$1" "$ff" | tr -d '\r')"
    local pycmd="print '""${str_to_search}""'.count('""$tag""')"
    local count=$(python -c "$pycmd")
    local tag_new="$2"
    if [ $count == 1 ] # if the substring was found only once in the line to be modified
    then
#	sed -i '' -e "${tag_line}"s/"${tag}"/"${tag_new}"/g "$ff" # Replace text at first line containing tag
	perl -i -pe "s/"$tag"/"${tag_new}"/" "$ff" # Replace text at first line containing tag
    else
	echo "Error: Value of" "$1" "'""$tag""'" "found more than once in" "$1" "line. Skipping modification of" "$1""."
    fi
}

# Print PHI fields before modification
echo "Old PHI:"
grep -a "Patient" "$ff" | grep -v "Position"
echo

# Modify PatientName and PatientID fields.
if [ ! -z "$nn" ] # if new PatientName was specified
then 
#    echo "Changing PatientName and PatientID to:" "$nn" 
    change_tag "PatientName" "$nn"
    change_tag "PatientID" "$nn"
else
#    echo "PatientName unchanged. PatientID changed to original PatientName."
    PatientName="$(grep -a "PatientName" "$ff" | cut -d " " -f 2- | cut -d $'\r' -f 1)"
    change_tag "PatientID" "$PatientName"
fi

# Modify remaining PHI fields
change_tag "PatientSex" "A"
change_tag "PatientBirthDate" "00000000"
change_tag "PatientAge" "999Y"
change_tag "PatientWeight" "0.000000"

# Print PHI fields before modification
echo "New PHI:"
grep -a "Patient" "$ff" | grep -v "Position"