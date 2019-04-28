#!/bin/bash

set -e
set -x

fslmaths aparc+aseg.nii.gz -uthr 2029 -thr 2029 aseg.visual.nii.gz 
flirt -in aseg.visual.nii.gz -ref $bold -out mask.nii.gz
