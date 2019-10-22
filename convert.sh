#!/bin/bash

set -x
set -e

echo "using $1 to convert freesurfer data"
mri_convert $1/mri/aparc+aseg.mgz aparc+aseg.nii.gz

echo "convert surfaces to vtk"
mkdir -p prf/surfaces
mris_convert --to-scanner $1/surf/lh.white prf/surfaces/lh.white.vtk
mris_convert --to-scanner $1/surf/rh.white prf/surfaces/rh.white.vtk
mris_convert --to-scanner $1/surf/lh.pial prf/surfaces/lh.pial.vtk
mris_convert --to-scanner $1/surf/rh.pial prf/surfaces/rh.pial.vtk
mris_convert --to-scanner $1/surf/lh.sphere prf/surfaces/lh.sphere.vtk
mris_convert --to-scanner $1/surf/rh.sphere prf/surfaces/rh.sphere.vtk
mris_convert --to-scanner $1/surf/lh.inflated prf/surfaces/lh.inflated.vtk
mris_convert --to-scanner $1/surf/rh.inflated prf/surfaces/rh.inflated.vtk
