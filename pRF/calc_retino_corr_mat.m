function calc_retino_corr_mat(templateDir,template,hemi,runs,func)

% Calculate the observed and predicted correlation matrices using a 
%   retinotopic template
%
%   Written by Andrew S Bock Dec 2015
%% set defaults
areas = {'V1' 'V2' 'V3'};

%% Get template and movie correlation matrices
ct = 0;
for ss = 1:length(areas)
    for tt = ss:length(areas)
        ct = ct + 1;
        areas_srctemplate = fullfile(templateDir,[hemi '.areas.' template '.nii.gz']);
        ecc_srctemplate = fullfile(templateDir,[hemi '.ecc.' template '.nii.gz']);
        pol_srctemplate = fullfile(templateDir,[hemi '.pol.' template '.nii.gz']);
        [corr_mat{ct},srcind,trgind] = ...
            make_model_corr(areas{ss},areas{tt},areas_srctemplate,ecc_srctemplate,pol_srctemplate);
        [movie_mat{ct}] = make_stimulus_corr(session_dir,runs,hemi,func,srcind,trgind);
    end
end