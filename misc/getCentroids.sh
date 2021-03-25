#!/bin/bash

atlas="$1"

# Get max value of atlas
max="$(fslstats "$atlas" -R | awk '{print $2}')"

# Print centre of mass for each label.
for ii in $(seq 1 $max)
do
    lower=$(python3 -c "print("$ii"-0.1)")
    upper=$(python3 -c "print("$ii"+0.1)")

    # Threshold out everything but current label, and print centre of gravity.
    cog="$(fslstats "$atlas" -l $lower -u $upper -C)"

    echo $ii $cog
done
