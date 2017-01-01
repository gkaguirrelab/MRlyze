function MP2RAGE_bkgnd(MP2RAGE,PD)
%   Remove noise MP2RAGE background, best if inputs are ACPC aligned
%
%   usage:
%   MP2RAGE_bkgnd(MP2RAGE,PD)
%
%   inputs:
%   MP2RAGE - path to MP2RAGE.nii.gz image (e.g.
%   <session_dir>/MP2RAGE/004/ACPC/MP2RAGE.nii.gz)
%   PD - path to proton-density-like image (e.g.
%   <session_dir>/MP2RAGE/002/ACPC/MP2RAGE.nii.gz)
%
%   output:
%   creates at file 'mask_50.nii.gz' in the PD directory
%   creates a file 'MP2RAGE_nobg.nii.gz' in the MP2RAGE directory
%
%   Written by Andrew S Bock Sept 2014

%% Find paths
[PDpath,~,~] = fileparts(PD);
[MP2RAGEpath,~,~] = fileparts(MP2RAGE);
%% Mask background in PD image 
% set any voxels with an intensity below 50 to 0
[~,~] = system(['fslmaths ' PD ' -thr 50 -bin ' ...
    fullfile(PDpath,'mask_50.nii.gz')]);
%% Apply mask to MP2RAGE
% Use the mask from the step above to remove the noise background
[~,~] = system(['fslmaths ' MP2RAGE ' -mul ' fullfile(PDpath,'mask_50.nii.gz') ...
    ' ' fullfile(MP2RAGEpath,'MP2RAGE_nobg.nii.gz')]);