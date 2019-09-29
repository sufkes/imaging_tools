## De-identify dicom files
#
# Usage:
# ./deid_dicom.sh <dicom file> <new PatientID and PatientName>
# 
# The PatientName and PatientID tags will be changed to the 
# string specified in the second argument.
#
# Requirements: dcmtk

file="$1";
name="$2";

# Change PatientName and PatientID to specified subject name.
# These tags will be changed if a new name was provided.
if [ ! -z "$name" ];
    then declare -a nametags=("(0010,0010)" "(0010,0020)")
    for nametag in "${nametags[@]}";
        do dcmodify -nb -q -i "${nametag}=${name}" "$file";
    done
fi