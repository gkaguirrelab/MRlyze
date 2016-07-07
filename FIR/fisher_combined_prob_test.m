function fisher_combined_prob_test(output_dir, session_dir,subject_name,subj_name, condition, runNums,func,thresh,inAnatomicalSpace, SUBJECTS_DIR)
% fisher_combined_prob_test(output_dir, session_dir,subject_name,subj_name, condition, runNums,func,thresh,inAnatomicalSpace, SUBJECTS_DIR)
% Loads pval.nii.gz or pval.anat.nii.gz from feat/stats directory and
% performs the Fisher's probability test.
%
%
% Input arguments:
% ================
%
%   output_dir : path to the output directory;
%   session_dir : path to the session directory;
%   subject_name : Freesurfer subject name;
%   subj_name : subject name;
%   condition : condition to include in the analysis;
%   runNums : number of bold runs to include in the analysis (must be in the same
%   session_dir);
%   func : smoothed functional data to include in the analysis (e.g.
%   wdrf.tf);
%   thresh : threshold for p-values
%   inAnatomicalSpace : if set to 'true' will apply Fisher's test in anatomical
%   space. Otherwise, Fisher's test will be applied in functional space 
%   SUBJECTS_DIR: Freesurfer subject Directory.
%
% Usage:
% ======
%
%   output_dir = '/data/jag/MELA/MelanopsinMR/Results';
%   session_dir = '/data/jag/MELA/MelanopsinMR/HERO_xxx1/060616'';
%   subject_name = 'HERO_xxx1_MaxMel';
%   subj_name = 'HERO_xxx1';
%   condition = 'MelPulses_400pct';
%   runNums = 1:9;
%   funcs = 'wdrf.tf';
%   inAnatomicalSpace = true;
%   fisher_combined_prob_test(output_dir, ...
%   session_dir,subject_name,subj_name, condition, ...
%   runNums,func,thresh,inAnatomicalSpace);
%
%
% 3/17/16   ms, gf      Written.
% 7/2/2016  gf, ms      Commented.

%
%
%  


%% set defaults
if ~exist('thresh','var')
    thresh = 0.05; % p-value threshold
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end

% Define the name of the volumes you want to apply Fisher's combined
% probability test to. This is likely in anatomical space, but we leave
% this up to the user here.
if inAnatomicalSpace
    pvalVolumeName = 'pval.anat.nii.gz';
else
    pvalVolumeName = 'pval.nii.gz';
end


%% get bold dirs
d = find_bold(session_dir);

%%
switch condition
    case {'MelPulses_400pct' , 'LMSPulses_400pct'}
% Get the feat stats dir
for i = runNums
    statsDir = fullfile(session_dir,d{i},[func '.feat'],'stats');
    dof = load(fullfile(statsDir,'dof'));
    pval = load_nifti(fullfile(statsDir,pvalVolumeName));
    
    tmp(:, :, :, i) = pval.vol;

end

% Take the natural log
tmp(tmp == 0) = NaN;
logTmp = log(tmp);
sumLogTmp = -2*sum(logTmp, 4);
pval.vol = sumLogTmp;
fprintf('\t * Saving out Fisher''s test...');
save_nifti(pval,fullfile(output_dir,[subj_name '_' condition '_' 'Fisher_Chisq.anat.nii.gz']));
fprintf('done!');

P = 1 - chi2cdf(sumLogTmp,2*length(runNums));
pval.vol = P;
fprintf('\t * Saving out Fisher''s test as p values...');
save_nifti(pval,fullfile(output_dir,[subj_name '_' condition '_' 'Fisher_pval.anat.nii.gz']));
fprintf('done!');

    case {'MaxMelCRF', 'MaxLMSCRF'}
         controls = {...
            '25pct'...
            '50pct'...
            '100pct'...
            '200pct'...
            '400pct'...
            'AttentionTask'...
            };