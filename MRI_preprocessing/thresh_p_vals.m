function thresh_p_vals(session_dir,runNums,func,thresh)

% Creates a binary mask of voxel where p-value < 'thresh' for ALL runs
% specified by 'runNums'. Also saves a 'ct' file, with values indicating
% the number of runs that pass threshold.
%
%   Usage:
%   thresh_p_vals(session_dir,runNums,func,thresh)
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
if ~exist('thresh','var')
    thresh = 0.05; % p-value threshold
end
%% get bold dirs
d = find_bold(session_dir);
%% Get the pvals from all the stats directories
ct = 0;
for i = runNums
    ct = ct + 1;
    statsDir = fullfile(session_dir,d{i},[func '.feat'],'stats');
    tmp = load_nifti(fullfile(statsDir,'pval.anat.nii.gz'));
    goodind(ct,:,:,:) = tmp.vol<thresh;
end
%% Create mask
mask = tmp; % create nifti structure
vals = tmp;
tmpmask = squeeze(sum(goodind));
vals.vol = tmpmask;
save_nifti(vals,fullfile(session_dir,[func '.p_vals_ct.nii.gz']));
maskind = tmpmask==length(runNums); % must be 1s for all runs
mask.vol = zeros(size(mask.vol));
mask.vol(maskind) = 1;
save_nifti(mask,fullfile(session_dir,[func '.p_vals_mask.nii.gz']));