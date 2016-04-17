function calc_cope_prop(session_dir,runNum,func)

% Calculates the proportion of cope1 and cope2
%
%   Usage:
%   calc_cope_prop(session_dir,runNums,func)
%
%   cope1.vol./(cope1.vol+cope2.vol)
%
%   Written by Andrew S Bock Nov 2015

%% Find bold directories
b = find_bold(session_dir);

%% Calc proportion
cope1 = load_nifti(fullfile(session_dir,b{runNum},[func '.feat'],'stats',...
    'cope1.anat.psc.nii.gz'));
cope2 = load_nifti(fullfile(session_dir,b{runNum},[func '.feat'],'stats',...
    'cope2.anat.psc.nii.gz'));
copesum = cope1;
copesum.vol = abs(cope1.vol)./(abs(cope1.vol)+abs(cope2.vol));
save_nifti(copesum,fullfile(session_dir,b{runNum},[func '.feat'],'stats',...
    'cope_1_2_prop.nii.gz'));
