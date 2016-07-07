function [means,sems] = psc_cope_get_means(session_dir,subject_name,runNums,func,ROIind,copeNames,startingCope)
% [means,sems] = psc_cope_get_means(session_dir,subject_name,runNums,func,ROIind,copeNames,starting Cope)
%
% Adapted from psc_cope
%
%
% Input arguments:
% ================
%
%   session_dir : path to the session directory;
%   subject_name : Freesurfer subject name;
%   runNums : number of bold runs to include in the analysis (must be in the same
%   func : smoothed functional data to include in the analysis (e.g.
%   wdrf.tf);
%   ROIind : index of ROI;
%   copeNames : names of the copes;
%   startingCope : from which cope the analysis starts (if not specified,
%   the analysis will start from the first cope).
%
% Usage:
% ======
%
%   session_dir = '/data/jag/MELA/MelanopsinMR/HERO_xxx1/060616'';
%   subject_name = 'HERO_xxx1_MaxMel';
%   runNums = 1:9;
%   func = 'wdrf.tf';
%   ROIind
%   copeNames = {...
%             'Sec00' ...
%             'Sec01' ...
%             'Sec02' ...
%             'Sec03' ...
%             'Sec04' ...
%             'Sec05' ...
%             'Sec06' ...
%             'Sec07' ...
%             'Sec08' ...
%             'Sec09' ...
%             'Sec10' ...
%             'Sec11' ...
%             'Sec12' ...
%             'Sec13' ...
%             };
%   startingCope = 15;
% [means,sems] = psc_cope_get_means(session_dir,subject_name,runNums,func,ROIind,copeNames,starting Cope)
%
%
% 4/5/16    ms, gf  Adapted.

%%
if ~exist('startingCope','var')
    startingCope = 1;
end
cc = startingCope:startingCope+length(copeNames)-1;

%% Find bold directories
b = find_bold(session_dir);

ct = 0;
for i = runNums
    ct = ct + 1;
    for j = 1:length(copeNames)
        fprintf('\n> Run %g: Calculating pct signal change for %s.\n', i, copeNames{j});
        copeout = fullfile(session_dir,b{i},[func '.feat'],'stats',['cope' num2str(cc(j)) '.anat.nii.gz']);
        meanout = fullfile(session_dir,b{i},[func '.feat'],'mean_func.anat.nii.gz');
        fprintf('\t * Loading %s...', copeout);
        ctmp = load_nifti(copeout);
        fprintf('done.\n');
        mtmp = load_nifti(meanout);
        fprintf('\t * Loading %s...', meanout);
        fprintf('done.\n');
        psctmp = (ctmp.vol./mtmp.vol)*100; % convert to percent signal change
        psctmp(psctmp==inf | psctmp==-inf) = nan; % set inf/-inf to nan
        eval([copeNames{j} '(ct) = nanmean(psctmp(ROIind));']);
        fprintf('> Pct signal change calculated.\n');
    end
end
for j = 1:length(copeNames)
    eval(['means(j) = nanmean(' copeNames{j} ');']);
    eval(['sems(j) = nanstd(' copeNames{j} ') / sqrt(length(' copeNames{j} '));']);
end