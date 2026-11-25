# pRFLife

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

**Fits population receptive field (pRF) models to BOLD time series and a stimulus image with mrTools, producing polar angle, eccentricity, receptive field width and variance-explained maps.**

This app will compute a population receptive field (i.e. pRF) analysis on a given time series. Inputs needed include: a time series nifti, a stimulus nifti, the measurements of the visual stimulus in space (stimimageunits), and dimensions of a mask (mask; optional). Outputs include: x.nii, y.nii, rfWidth.nii and r2.nii. Please see http://gru.stanford.edu/doku.php/mrTools/tutorialsprf for more detail.

Run it on brainlife.io: [https://doi.org/10.25663/bl.app.31](https://doi.org/10.25663/bl.app.31)

| Input | Datatype | Required | Description |
|---|---|---|---|
| `bold` | `neuro/func/task` | yes | BOLD fMRI time series from a retinotopic mapping experiment. |
| `stimimage` | `neuro/stimulus` | yes | Stimulus aperture image (NIfTI) shown during the BOLD run, as in the HCP retinotopy datasets. Should have one frame per BOLD volume. |
| `freesurfer` | `neuro/freesurfer` | yes | Freesurfer output used to create mask for visual cortex |

| Output | Datatype | Description |
|---|---|---|
| `prf` | `neuro/prf` | pRF maps (polarAngle, eccentricity, rfWidth and r2 NIfTI images) plus FreeSurfer white, pial, sphere and inflated surfaces as VTK for visualisation. |

| Parameter | Type | Default | Description |
|---|---|---|---|
| `visual_angle_width` | number | `16` | (in degree) full visual angle of the stimulation image. assumes central fixation. It should be 16 for HCP Retinotopy datasets. |
| `visual_angle_height` | number | `16` | (in degree) full visual angle of the stimulation image. assumes central fixation. It should be 16 for HCP Retinotopy datasets. |
| `frameperiod` | number | `1` | (in seconds) the time for each volume of the BOLD time series (typically the TR) |
| `prefitOnly` | boolean | `false` | Only do the grid search that finds approximately good parameters for the model and will not optimize with the nonlinear search. This will speed up the program but will not produce optimized fits. Nonetheless it should be good enough for debugging. |

Runtime: container `brainlife/freesurfer_on_mcr:6.0.0`; 4 CPUs; 120 GB memory; 24 h walltime; needs `FREESURFER_LICENSE`.

<!-- brainlife-app:end -->

![screen](screen.png)

This App compute pRF analysis from bold and HCP stimulus image.

# Author 

Justin Gardner




