## De-identify dicom files
#
# Usage:
# ./deid_dicom.sh <dicom file> [new PatientID and PatientName]
#
# If a single argument (the dicom file) is provided, the PatientName
# and PatientID tags will be unchanged (bad). If two arguments are 
# provided, the PatientName and PatientID tags will be changed to the 
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

## Replace PHI tags. These tags will be removed.
# AE tags - Application Entity; characters; 16 bytes maximum
declare -a tags_AE=("(0032,1021)" "(0040,0241)")
for tag in "${tags_AE[@]}";
    do dcmodify -nb -q -m "$tag=ANONYMIZED" "$file";
done

# AS tags - Age String; string with format nnnD, nnnW, nnnM, nnnY
declare -a tags_AS=("(0010,1010)")
for tag in "${tags_AS[@]}";
    do dcmodify -nb -q -m "$tag=000Y" "$file";
done

# CS tags - Code String; string with 16 bytes maximum
declare -a tags_CS=("(0010,0022)" "(0010,0035)" "(0010,0040)" "(0010,21A0)" "(0010,2298)")
for tag in "${tags_CS[@]}";
    do dcmodify -nb -q -m "$tag=ANONYMIZED" "$file";
done

# DA tags - Date; YYYYMMDD; 0-9, 8 bytes fixed
declare -a tags_DA=("(0010,0030)" "(0010,21D0)")
for tag in "${tags_DA[@]}";
    do dcmodify -nb -q -m "$tag=00000000" "$file";
done

# DS tags - Decimal String; floating point number (e.g. 0, 0.0, 1.03e+13) with 16 bytes maximum
declare -a tags_DS=("(0010,1020)" "(0010,1030)")
for tag in "${tags_DS[@]}";
    do dcmodify -nb -q -m "$tag=0" "$file";
done

# LO tags - Long String; characters with 64 characters maximum
declare -a tags_LO=("(0008,0080)" "(0008,1040)" "(0010,0021)" "(0010,0033)" "(0010,0034)" "(0010,1000)" "(0010,1040)" "(0010,1050)" "(0010,1080)" "(0010,1081)" "(0010,1090)" "(0010,2000)" "(0010,2110)" "(0010,2150)" "(0010,2152)" "(0010,21F0)" "(0010,2299)" "(0032,1020)" "(0032,1030)" "(0032,1033)" "(0038,0300)" "(0038,0400)" "(0038,0500)")
for tag in "${tags_LO[@]}";
    do dcmodify -nb -q -m "$tag=ANONYMIZED" "$file";
done

# LT tags - Long Text; characters with 10240 chars maximum
declare -a tags_LT=("(0008,4000)" "(0010,21B0)" "(0010,4000)" "(0032,4000)")
for tag in "${tags_LT[@]}";
    do dcmodify -nb -q -m "$tag=ANONYMIZED" "$file";
done

# PN tags - Person Name
declare -a tags_PN=("(0008,0090)" "(0008,009C)" "(0008,1048)" "(0008,1050)" "(0008,1060)" "(0008,1070)" "(0010,1001)" "(0010,1005)" "(0010,1060)" "(0010,2297)" "(0032,1032)" "(0040,0006)" "(0040,2008)" "(4008,010C)" "(4008,0114)" "(4008,0119)")
for tag in "${tags_PN[@]}";
    do dcmodify -nb -q -m "$tag=ANONYMIZED" "$file";
done

# SH tags - Short String
declare -a tags_SH=("(0008,0094)" "(0008,1010)" "(0010,2154)" "(0010,2160)" "(0010,2180)" "(0040,0010)" "(0040,0011)" "(0040,0242)" "(0040,0243)" "(0040,2009)" "(0040,2010)")
for tag in "${tags_SH[@]}";
    do dcmodify -nb -q -m "$tag=ANONYMIZED" "$file";
done

# SQ tags - Sequences; delete these.
declare -a tags_SQ=("(0008,0082)" "(0008,0096)" "(0008,009D)" "(0008,1049)" "(0008,1052)" "(0008,1062)" "(0008,1072)" "(0010,0050)" "(0010,0101)" "(0010,0102)" "(0010,1002)" "(0010,1021)" "(0032,1031)" "(0040,000B)")
for tag in "${tags_SQ[@]}";
    do dcmodify -nb -q -e "$tag=" "$file";
done

# ST tags - Short Text
declare -a tags_ST=("(0008,0081)" "(0008,0092)" "(0008,0116)")
for tag in "${tags_ST[@]}";
    do dcmodify -nb -q -m "$tag=ANONYMIZED" "$file";
done

# TM tags - Time HHMMSS.FFFFFF
declare -a tags_TM=("(0010,0032)")
for tag in "${tags_TM[@]}";
    do dcmodify -nb -q -m "$tag=000000.000" "$file";
done

# US tags - Unsigned Short; 2 bytes fixed.
declare -a tags_US=("(0010,21C0)")
for tag in "${tags_US[@]}";
    do dcmodify -nb -q -m "$tag=00" "$file";
done

#### Identifying tags to modify.
## Replace these tags with Study ID (e.g. SK040001)
#(0010,0010)    PatientName and 
#(0010,0020)    PatientID

## Replace these tags with anonymous values (e.g. 'ANONYMIZED', 00000000), or delete if type "Sequence" (SQ).
#(0008,0080)    Institution Name
#(0008,0081)    Institution Address
#(0008,0082)    Institution Code Sequence
#(0008,0090)    Referring Physician's Name
#(0008,0092)    Referring Physician's Address
#(0008,0094)    Referring Physician's Telephone Numbers
#(0008,0096)    Referring Physician Identification Sequence
#(0008,009C)    Consulting Physician's Name
#(0008,009D)    Consulting Physician Identification Sequence
#(0008,0116)    Responsible Organization
#(0008,1010)    Station Name
#(0008,1040)    Institutional Department Name
#(0008,1048)    Physician(s) of Record
#(0008,1049)    Physician(s) of Record Identification Sequence
#(0008,1050)    Performing Physician's Name
#(0008,1052)    Performing Physician Identification Sequence
#(0008,1060)    Name of Physician(s) Reading Study
#(0008,1062)    Physician(s) Reading Study Identification Sequence
#(0008,1070)    Operators' Name
#(0008,1072)    Operator Identification Sequence
#(0008,4000)    Identifying Comments (RET)
#(0010,0021)    Issuer of Patient ID
#(0010,0022)    Type of Patient ID
#(0010,0030)    Patient's Birth Date
#(0010,0032)    Patient's Birth Time
#(0010,0033)    Patient's Birth Date in Alternative Calendar
#(0010,0034)    Patient's Death Date in Alternative Calendar
#(0010,0035)    Patient's Alternative Calendar
#(0010,0040)    Patient's Sex
#(0010,0050)    Patient's Insurance Plan Code Sequence
#(0010,0101)    Patient's Primary Language Code Sequence
#(0010,0102)    Patient's Primary Language Code Modifier Sequence
#(0010,1000)    Other Patient IDs
#(0010,1001)    Other Patient Names
#(0010,1002)    Other Patient IDs Sequence
#(0010,1005)    Patient's Birth Name
#(0010,1010)    Patient's Age    
#(0010,1020)    Patient's Size    
#(0010,1021)    Patient's Size Code Sequence    
#(0010,1030)    Patient's Weight
#(0010,1040)    Patient's Address
#(0010,1050)    Insurance Plan Identification (RET)
#(0010,1060)    Patient's Mother's Birth Name
#(0010,1080)    Military Rank
#(0010,1081)    Branch of Service
#(0010,1090)    Medical Record Locator
#(0010,2000)    Medical Alerts
#(0010,2110)    Allergies
#(0010,2150)    Country of Residence
#(0010,2152)    Region of Residence
#(0010,2154)    Patient's Telephone Numbers
#(0010,2160)    Ethnic Group
#(0010,2180)    Occupation
#(0010,21A0)    Smoking Status
#(0010,21b0)    Additional Patient History
#(0010,21c0)    Pregnancy Status
#(0010,21d0)    Last Menstrual Date
#(0010,21f0)    Patient's Religious Preference
#(0010,2297)    Responsible Person
#(0010,2298)    Responsible Person Role
#(0010,2299)    Responsible Organization
#(0010,4000)    Patient Comments
#(0032,1020)    Scheduled Study Location (RET)
#(0032,1021)    Scheduled Study Location AE Title (RET)
#(0032,1030)    Reason for Study (RET)
#(0032,1031)    Requesting Physician Identification Sequence
#(0032,1032)    Requesting Physician
#(0032,1033)    Requesting Service
#(0032,4000)    Study Comments
#(0038,0300)    Current Patient Location
#(0038,0400)    Patient's Institution Residence
#(0038,0500)    Patient State
#(0040,0006)    Scheduled Performing Physician's Name
#(0040,000b)    Scheduled Performing Physician Identification Sequence
#(0040,0010)    Scheduled Station Name
#(0040,0011)    Scheduled Procedure Step Location
#(0040,0241)    Performed Station AE Title
#(0040,0242)    Performed Station Name
#(0040,0243)    Performed Location
#(0040,2008)    Order Entered By
#(0040,2009)    Order Enterer's Location
#(0040,2010)    Order Callback Phone Number
#(4008,010c)    Interpretation Author (RET)
#(4008,0114)    Physician Approving Interpretation (RET)
#(4008,0119)    Distribution Name (RET)
#(4008,011a)    Distribution Address (RET)
