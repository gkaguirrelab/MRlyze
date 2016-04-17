function mri_vol2vol(invol,targvol,outvol,reg,interp,inv)

%   Wrapper for Freesurfer's 'mri_vol2vol'
%
%   Usage:
%   mri_vol2vol(invol,targvol,outvol,reg,interp,inv)
%
% Resamples a volume into another field-of-view using various types of
% matrices (FreeSurfer, FSL, SPM, and MNI). This is meant to be used
% in conjunction with tkregister2.
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
if ~exist('interp','var')
    interp = 'nearest'; % could be 'trilin', 'nearest', or 'cubic'
end
if ~exist('inv','var')
    inv = 0; % inverse (default = off)
end
%% run mri_vol2vol
if inv
    system(['mri_vol2vol --mov ' invol ' --targ ' targvol ' --o ' outvol ...
        ' --reg ' reg ' --' interp ' --inv']);
else
    system(['mri_vol2vol --mov ' invol ' --targ ' targvol ' --o ' outvol ...
        ' --reg ' reg ' --' interp]);
end