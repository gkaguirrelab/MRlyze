function [means,sems] = psc_cope(session_dir,subject_name,runNums,func,ROIind,copeNames)

% Gets the mean and SEM percent signal change for an ROI across runs, for
%   each cope in a feat directory
%
%   Usage:
%   [means,sems] = psc_cope(session_dir,subject_name,runNums,func,ROIind,copeNames)
%
%   Written by Andrew S Bock Sep 2015

%% set defaults
if ~exist('func','var')
    func = 'wdrf.tf';
end
if ~exist('copeNames','var')
    copeNames = {...
        'Hz2' ...
        'Hz4' ...
        'Hz8' ...
        'Hz16' ...
        'Hz32' ...
        'Hz64' ...
        };
end    
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
%% Find bold directories
b = find_bold(session_dir);
%% Pull out beta values
ct = 0;
for i = runNums
    ct = ct + 1;
    for j = 1:length(copeNames)
        % Project cope to anatomical space
        copefile = fullfile(session_dir,b{i},[func '.feat'],'stats',['cope' num2str(j) '.nii.gz']);
        copeout = fullfile(session_dir,b{i},[func '.feat'],'stats',['cope' num2str(j) '.anat.nii.gz']);
        targ_vol = fullfile(SUBJECTS_DIR,subject_name,'mri','T1.mgz');
        tmpreg = listdir(fullfile(session_dir,b{i},'*bbreg.dat'),'files');
        bbreg_out_file = fullfile(session_dir,b{i},tmpreg{1}); % name registration file
        system(['mri_vol2vol --mov ' copefile ...
            ' --targ ' targ_vol ' --reg ' bbreg_out_file ...
            ' --o ' copeout ' --nearest']);
        % Project mean functional volume to anatomical space
        meanfile = fullfile(session_dir,b{i},[func '.feat'],'mean_func.nii.gz');
        meanout = fullfile(session_dir,b{i},[func '.feat'],'mean_func.anat.nii.gz');
        tmpreg = listdir(fullfile(session_dir,b{i},'*bbreg.dat'),'files');
        bbreg_out_file = fullfile(session_dir,b{i},tmpreg{1}); % name registration file
        system(['mri_vol2vol --mov ' meanfile ...
            ' --targ ' targ_vol ' --reg ' bbreg_out_file ...
            ' --o ' meanout ' --nearest']);
        % Calculate percent signal change (using copes)
        ctmp = load_nifti(copeout);
        mtmp = load_nifti(meanout);
        psctmp = (ctmp.vol./mtmp.vol)*100; % convert to percent signal change
        psctmp(psctmp==inf | psctmp==-inf) = nan; % set inf/-inf to nan
        eval([copeNames{j} '(ct) = nanmean(psctmp(ROIind));']);
    end
end
%% Create TTF
for j = 1:length(copeNames)
    eval(['means(j) = nanmean(' copeNames{j} ');']);
    eval(['sems(j) = nanstd(' copeNames{j} ') / sqrt(length(' copeNames{j} '));']);
end