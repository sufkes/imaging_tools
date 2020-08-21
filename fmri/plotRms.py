#!/usr/bin/env python

import os, sys
import argparse
from collections import OrderedDict
from matplotlib import pyplot as plt
import numpy as np

def plotRms(in_paths, out_dir, num_vols_to_exclude=20):
    for in_path in in_paths:
        #rel = in_path.replace("abs", "rel")
        rel = os.path.join(os.path.dirname(in_path), os.path.basename(in_path).replace("abs", "rel"))

        if (not os.path.exists(rel)):
            rel = None
            data_rel = None
        else:
            data_rel = np.loadtxt(rel)

            # Add a zero in the first slot so that volumes align
            data_rel_new = np.zeros(data_rel.shape[0] + 1 )
            data_rel_new[1:] = data_rel
            data_rel = data_rel_new


            
        # Find the best volumes to include if excluding exactly num_vol_to_exclude.
        # Volumes are considered worse if their RMS displacement relative to the previous volume is high.
        num_vols_to_include = len(data_rel) - num_vols_to_exclude # should be 100 for MagNUM

        ## Need to somehow assign a "badness" value to every volume (including the first and last volumes).
        # Badness method (1): Set the relative displacement of the first volume (volume 0) to half that of the second volume (volume 1). Set the relative displacement of the second volume to half it's actual value. The logic behind this is that th
        #data_rel_search = data_rel
        #data_rel_search[0] = data_rel_search[0]/2.0
        #data_rel_search[0] = data_rel_search[0]/2.0

        # Badness method (2): Set the rel_rms displacement of a volume to the average of the displacement from the last volume and the displacement to the next volume.
        data_rel_search = (data_rel + np.roll(data_rel, -1))/2.0
        
        best_mean_rel_rms = None
        for start_vol in range(num_vols_to_exclude+1):
            mean_rel_rms = data_rel_search[start_vol:start_vol+num_vols_to_include].mean()
            if (best_mean_rel_rms is None) or (mean_rel_rms < best_mean_rel_rms):
                best_start_vol = start_vol
                best_mean_rel_rms = mean_rel_rms
                
        first_vol = best_start_vol # volume number of first volume to include (first volume has volume number 0)
        last_vol = first_vol + num_vols_to_include - 1 # volume number of last volume to include

        # Extend the bounds in both directions as far as the mean relative displacement is less than 0.25 mm
        ext_first_vol = first_vol
        for test_vol in range(first_vol-1, -1, -1): # start at previous volume, go down to volume 0 (i.e. stop at volume -1), in steps of -1 volume.
            if data_rel_search[test_vol] <= 0.25:
                ext_first_vol = test_vol
                pass
            else:
                break # stop once we get to a volume with > 0.25 mm of relative RMS displacement.
        ext_last_vol = last_vol
        for test_vol in range(last_vol+1, len(data_rel_search), 1): # start at next volume, go up to last volume (i.e. stop at volume -1), in steps of -1 volume.
            if data_rel_search[test_vol] <= 0.25:
                ext_last_vol = test_vol
                pass
            else:
                break # stop once we get to a volume with > 0.25 mm of relative RMS displacement.

        # Print the number of volumes with rel_rms (average of motion values from previous and next frame).
        num_spikes = OrderedDict([(dist, 0) for dist in [0.25, 0.5, 1.0, 2.0, 4.0]])
        for vol in range(ext_first_vol, ext_last_vol+1):
            badness = data_rel_search[vol]
            for dist, count in num_spikes.iteritems():
                if (badness >= dist):
                    num_spikes[dist] += 1

        num_spikes_filename = "num_spikes.csv"
        num_spikes_path = os.path.join(out_dir, num_spikes_filename)

        # Write header if file does not exist.
        if (not os.path.exists(num_spikes_path)):
            with open(num_spikes_path, 'w') as handle:
                line = ",".join(["subject_id"] + [str(key)+"mm_spikes_in_included_vols" for key in num_spikes.keys()])+"\n"
                print line
                handle.write(line)
            
        with open(num_spikes_path, 'a') as handle:
            subject = os.path.basename(in_path).split('_')[0]
            line = ",".join([subject] + [str(val) for val in num_spikes.values()])+"\n"
            print line
            handle.write(line)
            
        data = np.loadtxt(in_path)

        fig = plt.figure()
        plt.plot(data)
        if (not data_rel is None):
            plt.plot(data_rel)
            
        plt.axis([0, 120, 0, 5])
        ref_vol = int(os.path.basename(in_path).split('_')[1][3:])
        plt.axvline(x=ref_vol, color='g') # line to indicate reference volume.
        plt.axvline(x=first_vol, color='r') # line to indicate start volume
        plt.axvline(x=last_vol, color='r') # line to indicate end volume
        plt.axvline(x=ext_first_vol, color='r', linestyle="--") # line to indicate start volume
        plt.axvline(x=ext_last_vol, color='r', linestyle="--") # line to indicate end volume
        title = "start: "+str(first_vol)+", end: "+str(last_vol)+"\nstart (ext): "+str(ext_first_vol)+", end (ext): "+str(ext_last_vol)
        plt.title(title)
        
        out_name = os.path.basename(in_path) + ".png"
        out_path = os.path.join(out_dir, out_name)
        fig.savefig(out_path)
        plt.close()

        # Print the important numbers to terminal.
        print os.path.basename(in_path).split("_")[0], ext_first_vol, ext_last_vol
        if (ref_vol < ext_first_vol) or (ref_vol > ext_last_vol):
            print "Warning: Reference volume falls outside the included range of volumes."
    

if (__name__ == '__main__'):
    # Create argument parser
    description = """Description of function"""
    parser = argparse.ArgumentParser(description=description)

    # Define positional arguments.
    parser.add_argument("in_paths", help="path to _abs.rms files generated from MCFLIRT, containing the RMS displacement from each volume to the reference volume. The script will automatically look for the _rel.rms files with analogous names, which are (I think) the displacement of each volume relative to the previous volume.", nargs="+")
    
    # Define optional arguments.
    parser.add_argument("-o", "--out_dir", help="directory to save figures in. Default is current directory", default=os.getcwd())
    
    # Parse arguments.
    args = parser.parse_args()

    # Do stuff
    plotRms(in_paths=args.in_paths, out_dir=args.out_dir)
