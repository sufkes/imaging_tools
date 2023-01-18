#!/bin/bash
##
usage () {
    echo 'Compute radial diffusivity image from eigenvalue images. Requires FSL function "fslmaths" to be available in the path.

Usage:
    getRadialDiffusivity.sh <path to second eigenvalue image (L2)> <path to third eigenvalue image (L3)> <path to output radial diffusivity image>'
}

test -z "$1" && { usage; exit 1; }

l2="$1"
l3="$2"
out_path="$3"
fslmaths "$l2" -add "$l3" -div 2 "$out_path"

				       
