function fisher_combined_prob_test(output_dir, session_dir,subject_name,subj_name, condition, runNums,func,thresh,inAnatomicalSpace, SUBJECTS_DIR)
% fisher_combined_prob_test(output_dir, session_dir,subject_name,subj_name, condition, runNums,func,thresh,inAnatomicalSpace, SUBJECTS_DIR)
%
% <GF> Please add input arguments and usage here
%
% Input arguments:
% ================
%
%   session_dir : 
%   ...
%
% Usage:
% ======
%
% <GF>
%
%
% 3/17/16   ms, gf      Written.
% 7/2/2016  gf, ms      Commented.

%
% Loads pval.nii.gz or pval.anat.nii.gz from feat/stats directory 
%
%   Usage:
%   fisher_combined_prob_test(session_dir,subject_name,runNums,func,thresh,SUBJECTS_DIR,inAnatomicalSpace)
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
%% Get the feat stats dir
for i = runNums
    statsDir = fullfile(session_dir,d{i},[func '.feat'],'stats');
    dof = load(fullfile(statsDir,'dof'));
    pval = load_nifti(fullfile(statsDir,pvalVolumeName));
    
    tmp(:, :, :, i) = pval.vol;

    
%     pval = fstat;
%     % Convert Fstats to p-values
%     pval.vol = 1-fcdf(fstat.vol,1,dof);
%     
%     
%     % Save volumes
%     save_nifti(pval,fullfile(statsDir,'pval.nii.gz'));
%     sigind = pval.vol<thresh;
%     pval.vol = zeros(size(pval.vol));
%     pval.vol(sigind) = 1;
%     save_nifti(pval,fullfile(statsDir,'pval.mask.nii.gz'));

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