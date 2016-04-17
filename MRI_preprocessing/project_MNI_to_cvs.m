function project_MNI_to_cvs(invol,outvol)

% Projects the FSL MNI vol (invol) to cvs_inMNI space (outvol)
%
%   Usage:
%   invol = '/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz';
%   outvol = '~/MNI1mm_in_cvs.nii.gz';
%   project_MNI_to_cvs(invol,outvol)
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
%% Register MNI to cvs
system(['tkregister2 --mov ' invol ' --targ ' ...
    fullfile(SUBJECTS_DIR,'cvs_avg35_inMNI152/mri/T1.mgz') ...
    ' --reg ' fullfile(SUBJECTS_DIR,'MNI2cvs.dat') ' --regheader --noedit']);
%% Apply registration (above)
system(['mri_vol2vol --mov ' invol ' --targ ' ...
    fullfile(SUBJECTS_DIR,'cvs_avg35_inMNI152/mri/T1.mgz') ...
    ' --reg ' fullfile(SUBJECTS_DIR,'MNI2cvs.dat') ' --o ' outvol]);