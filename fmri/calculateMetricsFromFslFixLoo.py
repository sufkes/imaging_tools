#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import pandas as pd

def getNegativeCases(cases_path):
    with open(cases_path, 'r') as f:
        noises = f.read()
    noises = noises.strip().lstrip('[').rstrip(']').replace(' ','').split(',')
    noises = [int(i) for i in noises]
    return noises

def getClasses(classes_path):
    negative_classes = []
    positive_classes = []
    with open(classes_path, 'r') as f:
        lines = f.readlines()[1:-1]
    for line in lines:
        if 'Unclassified Noise' in line:
            negative_classes.append(int(line.split(',')[0]))
        else:
            positive_classes.append(int(line.split(',')[0]))
    return negative_classes, positive_classes


def main(run_name, subjects, thresholds, in_dir, verbose=True):
    # Set up dataframe to store metrics for each threshold.
#    metrics = pd.DataFrame(columns=['threshold', 'TPR (mean)', 'TNR (mean)', '(3*TPR + TNR)/4 (mean)', 'accuracy (mean)'], index=thresholds)
    metrics = pd.DataFrame()
    metrics.index.name = 'threshold'
    for threshold in thresholds:
        tprs = []
        tnrs = []
        scores = [] # (3*TPR + TNR)/4, as calculated by fix.
        accs = []
        
        for subject in subjects:
            # Get path of the hand-labelled noise component file.
            in_name = 'hand_labels_noise.txt'
            subject_dir = os.path.join(in_dir, subject)
            cases_path = os.path.join(subject_dir, in_name)
            
            # Get path of the classifications file, e.g. MS040002_V01/fix4melview_training_LOO_new_V01_LOO_thr1.txt
            classes_path = os.path.join(subject, 'fix4melview_'+run_name+'_LOO_thr'+str(threshold)+'.txt')
            
            negative_cases = getNegativeCases(cases_path)
            negative_classes, positive_classes = getClasses(classes_path)
            ics = sorted(negative_classes + positive_classes)
            positive_cases = [i for i in ics if not i in negative_cases]

            true_positives = [i for i in positive_classes if i in positive_cases]
            false_positives = [i for i in positive_classes if i in negative_cases]
            true_negatives = [i for i in negative_classes if i in negative_cases]
            false_negatives = [i for i in negative_classes if i in positive_cases]
            
            # Numbers
            P = len(positive_cases)
            N = len(negative_cases)
            total = N + P

            TP = len(true_positives)
            FP = len(false_positives)
            TN = len(true_negatives)
            FN = len(false_negatives)
            accuracy = (TP + TN)/(P + N)
            
            TPR = TP/P
            FPR = FP/N
            TNR = TN/N
            FNR = FN/P

            tprs.append(TPR)
            tnrs.append(TNR)
            scores.append((3*TPR+TNR)/4)
            accs.append(accuracy)
            
            if verbose:
                print('-'*80)
                print(f'Subject: {subject}\nThreshold: {threshold}')
                print('')
                print(f'True positives: {true_positives}\nFalse positives: {false_positives}\nTrue negatives: {true_negatives}\nFalse negatives: {false_negatives}')
                print('')
                print(f'TPR = {TPR:.2f}\nTNR = {TNR:.2f}\nFPR = {FPR:.2f}\nFNR = {FNR:.2f}\naccuracy = {accuracy:.2f}')
#                print('')
#                print('TPR = TP/P = 1 - FNR\nTNR = TN/N = 1 - FPR\nFPR = FP/N = 1 - TNR\nFNR = FN/P = 1 - TPR\naccuracy = (TP+TN)/(P+N)')
#                print('')
#                print('Hand-labelled noise components (negative cases):')
#                with open(cases_path, 'r') as f:
#                    print(f.read())
#                print('')
#                print('Classifications file:')
#                with open(classes_path, 'r') as f:
#                    print(f.read())
        metrics.loc[threshold, 'TPR (mean)'] = np.mean(np.array(tprs))
        metrics.loc[threshold, 'TNR (mean)'] = np.mean(np.array(tnrs))
        metrics.loc[threshold, '(3*TPR + TNR)/4 (mean)'] = np.mean(np.array(scores))
        metrics.loc[threshold, 'accuracy (mean)'] = np.mean(np.array(accs))
        metrics.loc[threshold, 'TPR (median)'] = np.median(np.array(tprs))
        metrics.loc[threshold, 'TNR (median)'] = np.median(np.array(tnrs))
        metrics.loc[threshold, '(3*TPR + TNR)/4 (median)'] = np.median(np.array(scores))
        metrics.loc[threshold, 'accuracy (median)'] = np.median(np.array(accs))

    print('-'*80+'\n'+'-'*80)
    print('Leave-one-out summary statistics for each threshold:')
    print('')
    print(metrics.to_string())
    out_path = os.path.join(in_dir, 'metrics.csv')
    metrics.to_csv(out_path, index=True)
    ## Save summary LOO results.
    # e.g. ./custom_training_LOO_new_V01_LOO_results.txt
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """"""
    parser = argparse.ArgumentParser(description=description)
    
    # Define positional arguments.
    parser.add_argument('run_name', type=str, help='name of run')
    parser.add_argument('subjects', type=str, nargs='+', help='subject names')
    
    # Define optional arguments.
    parser.add_argument('-i', '--in_dir', default=os.getcwd(), help='run base directory')
    parser.add_argument('-t', '--thresholds', type=str, nargs='+', default=['1','2','5','10','20','30','40','50'], help='thresholds tested in leave one out (LOO) analysis')
    parser.add_argument('-v', '--verbose', action='store_true')

    # Print help if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))
