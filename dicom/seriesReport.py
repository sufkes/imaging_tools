#!/usr/bin/env python

import os, sys
import glob

in_dir = str(sys.argv[1]) # contains series folders with names like 001-series_description
check_up_to_num = 12

study_paths=glob.glob(in_dir+"/*")
study_paths.sort()

print "Listing studies which are missing series numbers <= "+str(check_up_to_num)
check_nums = range(1, check_up_to_num+1)
for study_path in study_paths:
    if study_path.endswith(".txt"):
        continue
    series_paths = glob.glob(study_path+"/*")
    series_paths.sort()
    series_nums = []
    for series_path in series_paths:
        series_num = int(os.path.basename(series_path).split("-")[0])
        series_nums.append(series_num)
    missing_nums = []
    for num in check_nums:
        if not num in series_nums:
            missing_nums.append(num)
    if len(missing_nums) > 0:
        print os.path.basename(study_path)+" missing series: "+", ".join([str(n) for n in missing_nums])
    
