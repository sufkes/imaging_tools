#!/bin/bash

ls ../../raw/ | cut -d . -f 1  | while read ss;
do
        for ii in $(seq 0 119);
        do
                mean=after_mcflirt_allvols/"$ss"_ref"$ii"_abs_mean.rms;
                all=after_mcflirt_allvols/"$ss"_ref"$ii"_abs.rms;
                echo -n "$ss" "$ii" "$(cat "$mean")" "$(fmriMotionReport.py "$all" | grep Status | awk '{print $2}')"; echo;
        done;
done
