function save_psc_cope(session_dir,ROI,hemi,chiSqVol,featDirs,copeNums,outFile,threshVal)

% Runs 'psc_cope' and saves the outputs to a file
%
%   Usage:
%   save_psc_cope(outFile,chiSqVol,roiVol,featDirs,copeNums)
%
%   Written by Andrew S Bock May 2016

%% Set defaults
if ~exist('threshVal','var')
    threshVal = 0.05;
end
%% Load inputs, threshold using chi-squared p-values and ROI mask
% Get areas and ecc
areas = load_nifti(fullfile(session_dir,'anat_templates',...
    [hemi '.areas.anat.vol.nii.gz']));
ecc = load_nifti(fullfile(session_dir,'anat_templates',...
    [hemi '.ecc.anat.vol.nii.gz']));
lgn = load_nifti(fullfile(session_dir,'anat_templates',[hemi '.LGN.nii.gz']));
sc = load_nifti(fullfile(session_dir,'anat_templates',[hemi '.SC.nii.gz']));
chiSquared = load_nifti(chiSqVol);
switch ROI
    case 'V1all'
        ROIind = find(abs(areas.vol)==1 & chiSquared.vol<threshVal); % all of V1
    case 'V1low'
        ROIind = find(abs(areas.vol)==1 & ecc.vol<=5 & chiSquared.vol<threshVal);
    case 'V1mid'
        ROIind = find(abs(areas.vol)==1 & (ecc.vol>5 & ecc.vol<=15) & chiSquared.vol<threshVal);
    case 'V1high'
        ROIind = find(abs(areas.vol)==1 & (ecc.vol>15 & ecc.vol<=40) & chiSquared.vol<threshVal);
    case 'V2andV3'
        ROIind = find((abs(areas.vol)==2 | abs(areas.vol)==3) & chiSquared.vol<threshVal);
    case 'LGN'
        ROIind = find(lgn.vol & chiSquared.vol<threshVal);
    case 'SC'
        ROIind = find(sc.vol & chiSquared.vol<threshVal);
end
%% Run 'psc_cope'
[means,stds,sems] = psc_cope(featDirs,ROIind,copeNums);

%% Save output
save(outFile,'means','stds','sems','session_dir','ROI','hemi','chiSqVol','featDirs','copeNums','ROIind');