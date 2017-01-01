function [meanStats] = mean_surface_stats(session_dir,hemi,statName,numStats)

% Calculates the mean statistic across runs on the cortical surface
%
%   Usage:
%   [meanStats] = mean_surface_stats(session_dir,hemi,statName,numStats)
%
%   Written by Andrew S Bock Sep 2015

%% set defaults
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
%% Find bold directories
b = find_bold(session_dir);
%% Pull out beta values
for i = 1:numStats
    for j = 1:length(b)
        statFile = fullfile(session_dir,b{j},'sdbrf.tf.feat','stats',[statName num2str(i) '.nii.gz']);
        outname = fullfile(session_dir,b{j},'sdbrf.tf.feat','stats',[statName num2str(i) '.surf.' hemi '.nii.gz']);
        if ~exist(outname,'file');
            bbreg_out_file = fullfile(session_dir,b{j},'brf_bbreg.dat'); % name registration file
            [~,~] = system(['mri_vol2surf --src ' statFile ...
                ' --reg ' bbreg_out_file ' --hemi ' hemi ...
                ' --out ' outname ' --projfrac 0.5']);
        end
        tmp = load_nifti(outname);
        stats(i).vol(:,j) = tmp.vol;
    end
end
%% Calculate mean for each stat
for i = 1:numStats
    meanStats(:,i) = mean(stats(i).vol,2);
end