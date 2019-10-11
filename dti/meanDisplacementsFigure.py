#!/usr/bin/env python

import os, sys
import argparse

from matplotlib import pyplot as plt
import numpy as np

from pprint import pprint
from glob import glob

def meanDisplacementFigure(in_dir, out_path, plot_title, map_path):
    def getIDs(warp_paths):
        ids = set()
        for path in warp_paths:
            id = path.split('/')[-1].split('_to_')[0] # gets the ID of subject whose image is being registered.
            ids.add(id)
        ids = sorted(list(ids))
        return ids

    def makeScoreMatrix(warp_paths, ids):
        num_ids = len(ids)
        score_matrix = np.zeros(shape=(num_ids,num_ids))

        ## Add scores to matrix.
        for warp_index in range(len(warp_paths)):
            warp_path = warp_paths[warp_index]
            with open(warp_path, 'rb') as handle:
                lines = handle.readlines()
            mean_string = list(lines)[0].split()[0]
            mean = np.float(mean_string)

            # Figure out the indices of the matrix element to which the current score corresponds (should be in order anyway).
            id_from = warp_path.split('/')[-1].split('_to_')[0]
            id_to = warp_path.split('/')[-1].split('_to_')[1].split('_warp')[0]
            index_from = ids.index(id_from)
            index_to = ids.index(id_to)

            # Add current score to matrix.
            score_matrix[index_from, index_to] = mean

            sys.stdout.write('\r')
            sys.stdout.write('%.2f%% complete' % (float(warp_index)/float(len(warp_paths))*1e2,))
            sys.stdout.flush()

        sys.stdout.write('\r')
        sys.stdout.write('%.2f%% complete' % (float(100),))
        sys.stdout.write('\n')
        sys.stdout.flush()

        return score_matrix

    def makeFigure(score_matrix, out_path, plot_title):
        plt.figure()
        plt.imshow(score_matrix, origin='upper', vmin=0, vmax=20, cmap='jet')
        plt.colorbar()
        plt.title(plot_title,y=1.02)
        plt.xlabel("Subject image as target")
        plt.ylabel("Subject image for registration")

        plt.tight_layout()
        
        plt.savefig(out_path, bbox_inches='tight')
        plt.close()       
        return

    def makeMap(map_path, ids):
        if (map_path == None):
            return
        with open(map_path, 'wb') as handle:
            for subject_number in range(len(ids)):
                subject_id = ids[subject_number]
                line = str(subject_number) + ', ' + subject_id
                handle.write(line+'\n')

    # Find all of the warp files and generate a list of subjects
    warp_paths = sorted(glob(in_dir+"/*_to_*warp.msf"))
    ids = getIDs(warp_paths)
    
    # Make matrix that stores all of the mean displacement scores
    score_matrix = makeScoreMatrix(warp_paths, ids)
    
    # Create the color plot and ID map file.
    makeFigure(score_matrix, out_path, plot_title)
    makeMap(map_path, ids)

    return

if (__name__ == '__main__'):
    # Create argument parser
    description = """Generate Figure 1 (A & B) from paper: G. Ball et al. An optimised tract-based spatial statistics protocol for neonates: Applications to prematurity and chronic lung disease. NeuroImage 53 (2010).

Figure shows mean displacements for subject-to-target registration steps in TBSS."""
    epilog = None # """Text to follow argument explantion """
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("in_dir", help="path to TBSS FA directory containing of the form SubjectA_FA_to_SubjectB_FA_warp.msf, which indicate the mean (and median) displacement for each registration.", type=str)
    
    # Define optional arguments.
    parser.add_argument('-o', '--out_path', help="path to directory where figure should be saved. Default: current working directory", type=str, default=os.path.join(os.getcwd(), 'plot.png'))
    parser.add_argument('-t', '--plot_title', help='Title of plot', type=str, default='')
    parser.add_argument('-m', '--map_path', help='Path to save number-to-Subject ID map to. By default, do not save the map.', type=str)

    # Print help if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Do stuff.
    meanDisplacementFigure(args.in_dir, args.out_path, args.plot_title, args.map_path)
