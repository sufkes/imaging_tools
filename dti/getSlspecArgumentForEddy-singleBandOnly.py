#!/usr/bin/env python3

import os
import sys
import argparse
import json

def main(json_path, out_path, num_slices, save_guesses):
    # Open BIDS sidecar JSON file.
    with open(json_path, 'r') as handle:
        sidecar = json.load(handle)

    try:
        slice_timing = sidecar['SliceTiming'] # Read the slice timings.
        slice_timing_indexed = list(enumerate(slice_timing)) # Index each slice time (zero-based).
        slice_timing_indexed.sort(key=lambda x:x[1]) # Sort based on time.
        slice_orders = [pair[0] for pair in slice_timing_indexed] # Grab just the orders, removing the times.
        slice_timing_known = True
    except KeyError:
        slice_timing_known = False
        print("Warning: SliceTiming tag not found in BIDS sidecar. To save guesses for the slice ordering, run with flag --save_guesses.")

    if (num_slices is None):
        try:
            num_slices = len(slice_timing)
        except:
            raise Exception('Could not determine number of slices in image; rerun with num_slices set.')
    
        
    # Make two guesses of the slice ordering, assuming interleaved acquisition with the 
    slice_orders_guess_0 = [None] * num_slices
    slice_orders_guess_1 = [None] * num_slices

    stop0 = (num_slices - 1) // 2 + 1
    slice_orders_guess_0[0:stop0] = range(0,num_slices, 2)
    slice_orders_guess_0[stop0:] = range(1, num_slices, 2)

    stop1 = (num_slices) // 2
    slice_orders_guess_1[0:stop1] = range(1,num_slices, 2)
    slice_orders_guess_1[stop1:] = range(0, num_slices, 2)

    if num_slices % 2 == 0:
        slice_orders_guess_final = slice_orders_guess_1
    else:
        slice_orders_guess_final = slice_orders_guess_0
    
    if slice_timing_known:
        # Check if I can guess exactly right for SK/MS scans.
        guess_correct = (slice_orders == slice_orders_guess_final)
        if not guess_correct:
            print('Incorrect guess for file:', json_path)
            print('Number of slices:', num_slices)
            print('Actual:')
            print(slice_orders)
            print('Guess:')
            print(slice_orders_guess_final)
            print('')
        
        # Write to file.
        with open(out_path, 'w') as handle:
            handle.writelines([str(slice_num) + '\n' for slice_num in slice_orders])

    if save_guesses:
        out_dir= os.path.dirname(out_path)
        guess_out_path_0 = os.path.join(out_dir, 'my_slspec-guess0.txt')
        guess_out_path_1 = os.path.join(out_dir, 'my_slspec-guess1.txt')
        guess_out_path_final = os.path.join(out_dir, 'my_slspec-guess_final.txt')
        with open(guess_out_path_0, 'w') as handle:
            handle.writelines([str(slice_num) + '\n' for slice_num in slice_orders_guess_0])
        with open(guess_out_path_1, 'w') as handle:
            handle.writelines([str(slice_num) + '\n' for slice_num in slice_orders_guess_1])
        with open(guess_out_path_final, 'w') as handle:
            handle.writelines([str(slice_num) + '\n' for slice_num in slice_orders_guess_final])
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Generate the text file to be input as the --slspec argument to FSL eddy. It is a text file denoting the temporal order in which slices were acquired.
NB: This only handles the single-channel acquisitions."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("json_path", type=str, help="BIDS sidecar JSON file generated by dcm2niix. Must contain the SliceTiming array.")
    parser.add_argument("out_path", type=str, help="Path of output file which contains slice ordering to be used as the --slspec argument in FSL eddy.")
    
    # Define optional arguments.
    parser.add_argument("-n", "--num_slices", type=int, help="number of slices in corresponding image. This will be used to make guesses of the slice ordering, assuming interleaved acquisition, then either compare those guesses with the actual ordering, or save the guess in the case the the actual ordering is unknown due to a missing SliceTiming entry in the JSON file.")
    parser.add_argument('-g', '--save_guesses', action='store_true', help='guess the slice timing and save the guesses to file.')
    
    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))
