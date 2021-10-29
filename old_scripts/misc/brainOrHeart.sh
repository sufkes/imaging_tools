#!/bin/bash

in_dir="$1"


## Define brain keywords
brain_string='sss|dti|diffusion|ax-t1-3d|cor-t1-3d|sag-t1-3d|ax-3d-t1|cor-3d-t1|sag-3d-t1|ax-t2|cor-t2|sag-t2|ax-t1-mpr|cor-t1-mpr|sag-t1-mpr|ax-t1-tse|cor-t1-tse|sag-t1-tse|ax-t1-flash|cor-t1-flash|sag-t1-flash|ax-3d-mpr|cor-3d-mpr|sag-3d-mpr|ax-ssfp-3d|cor-ssfp-3d|sag-ssfp-3d|ax-3d-ssfp|cor-3d-ssfp|sag-3d-ssfp|svs|mag_image|pha_image|mip_image|swi_image'
# superior saggital sinus (SSS)

## Define heart keywords
heart_string='ventricles|svc|mpa|lpa|rpa|dao|aao|uv|2ch|3ch|4ch' # tf2d15_2ch, tf2d15_4ch pc, tof, carotid, neck, vessel? -- consider neck to be brain?
# superior vena cava (SVC)
# main, right and left pulmonary arteries (MPA, RPA and LPA)
# descending aorta (DAO)
# ascending aorta (AAO)
# umbilical vein (UV)? 
# 2-chamber (2CH)
# 3-chanber (3CH)
# 4-chamber (4CH)

## Neutral terms.
# integrated parallel imaging techniques (iPAT) -- iPat is the term Siemens uses for its parallel imaging implementation. It stands for integrated parallel imaging techniques.
# short-axis (SAX)
# Magnetic Resonance Venography (MRV)
# Magnetic Resonance Angiography (MRA)

# 'Standard' cardiac views for the left ventricle (LV) typically comprise 2-chamber, 3-chamber, 4-chamber, short-axis (SAX) views, or coronal LVOT. Right ventricle (RV) views typically comprise RVOT and RV 2-chamber views.

# Brain series
#echo "Brain series:"
#ls "$in_dir" | grep -Ei "$brain_string" | grep -Eiv "$heart_string"
num_brain=$(ls "$in_dir" | grep -Ei "$brain_string" | grep -Eiv "$heart_string" | wc -l)

# Heart series
#echo "Heart series:"
#ls "$in_dir" | grep -Ei "$heart_string" | grep -Eiv "$brain_string"
num_heart=$(ls "$in_dir" | grep -Ei "$heart_string" | grep -Eiv "$brain_string" | wc -l)

# Unknown
#echo "Unknown series:"
#ls "$in_dir" | grep -Eiv "$heart_string" | grep -Eiv "$brain_string"

if [ $num_brain -gt 0 ]
then
    if [ $num_heart -gt 0 ]
    then
	scan_type="C"
    else
	scan_type="B"
    fi
elif [ $num_heart -gt 0 ]
then
    scan_type="H"
else
    scan_type="N"
fi

#echo "Scan type:" $scan_type
echo $scan_type
    
    
	   
	      
