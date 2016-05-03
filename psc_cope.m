function [means,stds,sems] = psc_cope(featDirs,ROIind,copeNums)

% Gets the mean, std. deviation, and SEM percent signal change for an ROI 
%   across runs
%
%   Usage:
%   [means,stds,sems] = psc_cope(featDirs,ROIind,copeNums)
%
%   Inputs:
%   featDirs = cell of paths to FEAT directories
%   ROIind = index into a region of interest (anatomical space)
%   copeNums = vector of cope numbers 
%   (e.g. copeNums = 1:6; % cope1.anat.nii.gz ... cope6.anat.nii.gz)
%
%   Written by Andrew S Bock Sep 2015

%% Compute percent signal change
pscCope = cell(1,length(copeNums));
for i = 1:length(featDirs)
    % mean functional volume in anatomical space
    meanout = fullfile(featDirs{i},'mean_func.anat.nii.gz');
    for j = 1:length(copeNums)
        % copes in anatomical space
        copeout = fullfile(featDirs{i},'stats',['cope' num2str(copeNums(j)) '.anat.nii.gz']);
        % Calculate percent signal change (using copes)
        ctmp = load_nifti(copeout);
        mtmp = load_nifti(meanout);
        psctmp = (ctmp.vol./mtmp.vol)*100; % convert to percent signal change
        psctmp(psctmp==inf | psctmp==-inf) = nan; % set inf/-inf to nan
        pscCope{j}(i) = nanmean(psctmp(ROIind));
    end
end
%% Create TTF
means = nan(1,length(copeNums));
stds = nan(1,length(copeNums));
sems = nan(1,length(copeNums));
for j = 1:length(copeNums)
    means(j) = nanmean(pscCope{j});
    stds(j) = nanstd(pscCope{j});
    sems(j) = nanstd(pscCope{j}) / sqrt(length(pscCope{j}));
end