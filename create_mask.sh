#!/bin/bash

set -e
set -x

echo "using $1"

fslmaths aparc+aseg.nii.gz -uthr 2029 -thr 2029 aseg.visual.nii.gz 
flirt -in aseg.visual.nii.gz -ref $1 -out mask.nii.gz
