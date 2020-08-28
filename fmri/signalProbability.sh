#!/usr/bin/env bash

include_prob_map=include_2mm.nii.gz

find . -mindepth 3 -mindepth 3 -name *probmap* | sort | while read ff
do
    # Get mean probability (of nonzero voxels) from ICA probmap. This is used as a proxy for the "volume" of the ICA.
    mean_prob_map="$(fslstats "$ff" -M)"
    #echo $mean_prob_map

    # Get an image of the (probably of being in a brain region of interest -- cortex, DGM, brainstem) * (probabilty from ICA probmap)
    temp="$(mktemp)"
    fslmaths "$ff" -mul "$include_prob_map" "$temp"
    mean_prob2="$(fslstats "$temp" -M)"
    #echo $mean_prob2
    
    # Get the measure for the probability that the component is signal.
    val=$(python -c "print ${mean_prob2}/${mean_prob_map}")
    echo $val > "$(dirname "$ff")"/../signal_prob"$(basename "$ff" | cut -d _ -f 2 | cut -d "." -f 1)".txt
    rm "$temp"
done
