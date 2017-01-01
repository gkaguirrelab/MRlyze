function apply_B0(mag,fieldmap,epi,echospacing,warpdir,smoothing,TR)
%   Applies the B0 field map to epi data. Assumes the epi data has already
%   been motion corrected. If it has not, use my function "mcflirt.m"
%
%   Usage;
%   apply_B0(mag,fieldmap,epi,echospacing,smoothing)
%
%   defaults:
%   mag = no default, must specify (do not include .nii.gz
%   fieldmap = no default, must specify
%   epi = no default, must specify
%   echospacing = 0.64506; % from Siemens EPI at Penn, aka Dwell Time
%   warpdir = y; warp direction is "y"
%   smoothing = 0; % don't smooth fieldmap, other values are FWHM (mm)
%   TR = 0; choose the epi TR to use for alignment with the B0 images
%   (alternatively could use the middle TR, or any other TR)
%
%   Written by Andrew S Bock Sept 2014

%% Set default parameters
if ~exist('mag','var')
    error('no "mag" image defined');% must define a magnitude image
end
if ~exist('fieldmap','var')
    error('no "fieldmap" image defined');% must define a fieldmap image
end
if ~exist('epi','var')
    error('no "epi" image defined');% must define an epi volume
end
if ~exist('echospacing','var')
    echospacing = 0.64506; % Echo Spacing (sec) (from HUP6 at UPenn)
end
if ~exist('warpdir','var')
    warpdir = y; % warp direction (from HUP6 at UPenn)
end
if ~exist('smoothing','var')
    smoothing = 0; % don't do smoothing
end
if ~exist('TR','var')
    TR = 0; % use the 1st TR image for alignment with the B0 images
end
%% Brain extract magnitude and fieldmap images
[~,magname,~] = fileparts(mag);
[~,fieldmapname,~] = fileparts(fieldmap);
system(['bet ' mag ' ' magname '_brain -f .3 -m']);
system(['fslmaths '  fieldmap ' -mas ' magname '_mask ' fieldmapname '_brain.nii.gz']);
%% Smoothing
if smoothing
    system(['fugue --loadfmap=' fieldmapname '_brain -s ' smoothing ...
        ' --savefmap=' fieldmapname '_brain.nii.gz']);
end
%% Extract middle time point from EPI data
[~,epiname,~] = fileparts(epi);
system(['fslroi ' epiname ' ' epiname '_1TR ' TR ' 1']);
system(['bet ' epiname '_1TR ' epiname '_brain -m -f .3']);
fugue -v -i FieldMap/Magnitude_brain --unwarpdir=y --dwell=0.000720 --loadfmap=FieldMap/FieldMap.nii.gz -w FieldMap/Magnitude_brain_warpped
flirt -in FieldMap/Magnitude_brain_warpped.nii.gz -ref DTI/nodif_brain.nii.gz -out FieldMap/Magnitude_brain_warpped_2_nodif_brain -omat FieldMap/fieldmap2diff.mat
flirt -in FieldMap/FieldMap_brain_s4.nii.gz -ref DTI/nodif_brain -applyxfm -init FieldMap/fieldmap2diff.mat -out FieldMap/FieldMap_brain_s4_2_nodif_brain
fugue -v -i DTI/data1_corr.nii.gz --icorr --unwarpdir=y --dwell=0.000720 --loadfmap=FieldMap/FieldMap_brain_s4_2_nodif_brain.nii.gz -u data/data.nii.gz
