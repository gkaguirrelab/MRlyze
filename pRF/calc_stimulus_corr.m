function [template_corr] = calc_stimulus_corr(session_dir,runs,hemi,template,src_area,trg_area,cluster)

% Calculates the correlation between stimulus cross-correlation matrix, and
% a predicted cross-correlation matrix
%
%   Usage:
%   [template_corr] = calc_stimulus_corr(session_dir,runs,hemi,template,src_area,trg_area,cluster)
%   Written by Andrew S Bock Aug 2015

%% Set defaults
if ~exist('cluster','var')
    cluster = 0; % assume not running on cluster
end
%% Set pRF_dir
pRF_dir = fullfile(session_dir,'pRFs','model_templates','decimated_templates');

%% Set template file names
areas_srctemplate = fullfile(pRF_dir,[hemi '.areas_' template '.nii.gz']);
ecc_srctemplate = fullfile(pRF_dir,[hemi '.ecc_' template '.nii.gz']);
pol_srctemplate = fullfile(pRF_dir,[hemi '.pol_' template '.nii.gz']);
[corr_mat,srcind,trgind] = ...
    make_model_corr(src_area,trg_area,areas_srctemplate,ecc_srctemplate,pol_srctemplate);
[movie_mat] = make_stimulus_corr(session_dir,runs,hemi,srcind,trgind);
template_corr = corr(corr_mat(:),movie_mat(:));
%% If running on cluster, save output value to text file
if cluster
    [~,~] = system(['echo ' num2str(template_corr) ' > ' fullfile(pRF_dir,...
        [hemi '.' template '.corr.txt'])]);
end