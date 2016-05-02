function computeFisherChisq(inVols,outDir,outName)
% Uses Fisher's method to combine p-values from several input volumes into
% one test statistic.
%
%   Usage:
%   computeFisherChisq(inVols,outDir,outName)
%
%   Example:
%   % Get the input volumes
%   for i = 1:10
%       inVol{i} = '/some/directory/name{i}/pval.anat.nii.gz';
%   end
%   outDir = '/some/out/directory';
%   outName = 'subjectName_conditionName_funcName';
%   computeFisherChisq(inVols,outDir,outName);

%% Load in the input p-value volumes
for i = 1:length(inVols)
    pval = load_nifti(inVols{i});
    tmp(:, :, :, i) = pval.vol;
end
%% Take the natural log
tmp(tmp == 0) = NaN;
logTmp = log(tmp);
sumLogTmp = -2*sum(logTmp, 4);
%% Save the Chi-squared output
out = tmp;
out.vol = nan(size(out.vol));
out.vol = sumLogTmp;
disp(['Saving ' fullfile(outDir,[outName '_Fisher_Chisq.anat.nii.gz'])]);
save_nifti(out,fullfile(outDir,[outName '_Fisher_Chisq.anat.nii.gz']));
disp('done!');
%% Save the p-value output
out = tmp;
out.vol = nan(size(out.vol));
out.vol = 1 - chi2cdf(sumLogTmp,length(inVols));
disp(['Saving ' fullfile(outDir,[outName '_Fisher_pval.anat.nii.gz'])]);
save_nifti(out,fullfile(outDir,[outName '_Fisher_pval.anat.nii.gz']));
disp('done!');