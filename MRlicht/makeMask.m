function makeMask(inVol,outVol,thresh)

% Saves an output volume that is a binary mask of the input volume, using a
% quantile threshold
%
%   Usage:
%   makeMask(inVol,outVol,thresh)
%
%   Example:
%   inVol = /path/to/my/stat/volume.nii.gz;
%   outvol = /path/to/my/stat/mask.nii.gz;
%   thresh = 0.9;
%   makeMask(inVol,outVol,thresh);
%
%   Written by Andrew S Bock Mar 2016

%% Create the mask
% Load input volume
In = load_nifti(inVol);
% Threshold
Y = quantile(In.vol(:),thresh);
% Find voxels that pass threshold
tInd = In.vol>Y;
% Save Mask
tMask = In; 
tMask.vol = zeros(size(tMask.vol));
tMask.vol(tInd) = 1;
save_nifti(tMask,outVol);