function plot_fine_template_varexp(session_dir,subject_name,hemi,template,makeplot)

% Plots the output of 'calc_fine_template_varexp' on the cortical surface
%
%   Usage:
%   plot_template_varexp(session_dir,subject_name,hemi,template,makeplot)
%
%   Written by Andrew S Bock Oct 2015

%% Set defaults
pRF_dir = fullfile(session_dir,'pRFs','fine_model_templates','decimated_templates');
names = {'V1' 'V2' 'V3' 'all'};
if ~exist('makeplot','var')
    makeplot = 0; % plot using surface_plot
end
%% Find areas
areas = load_nifti(fullfile(pRF_dir,[hemi '.areas.' template '.nii.gz']));
V1ind = find(areas.vol == -1 | areas.vol == 1);
V2ind = find(areas.vol == -2 | areas.vol == 2);
V3ind = find(areas.vol == -3 | areas.vol == 3);
%% Load in error .mat file
load(fullfile(pRF_dir,[hemi '.' template '.varexp.mat']));
%% Plot on surface
src_surf = '0.1.inflated';
trg_surf = 'inflated';
tmp = areas;
tmp.vol = nan(size(tmp.vol));
% The 4th row is the entire matrix (1=V1,2=V2,3=V3)
for i = 1:4
    tmp.vol(V1ind) = var_exp(i,V1inds);
    tmp.vol(V2ind) = var_exp(i,V2inds);
    tmp.vol(V3ind) = var_exp(i,V3inds);
    inFile = fullfile(pRF_dir,[hemi '.' template '.varexp.' names{i} '.0.1.inflated.nii.gz']);
    outFile = fullfile(pRF_dir,[hemi '.' template '.varexp.' names{i} '.inflated.nii.gz']);
    save_nifti(tmp,inFile);
    decimate_surf(subject_name,hemi,src_surf,trg_surf,inFile,outFile);
    if makeplot
        surface_plot('var',outFile,subject_name,hemi);
    end
end