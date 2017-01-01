function FSL2Freesurfer(input_vol,out_path,outname,FWHM)

% Takes in an input volume, outputs that volume on the cortical surface
%   Assumptions:
%       input volume - in FSL's 2mm MNI space
%       output surface - Freesurfer's cvs_avg35_inMNI152 space
%
%   Usage:
%       FSL2Freesurfer(input_vol,out_path,outname,FWHM)
%
%   Example:
%       FSL2Freesurfer('~/zstat1.nii.gz','~/data/','zstat1_surf',5);
%
%   Written by Andrew S Bock Jun 2015

%% Set defaults
if ~exist('FWHM','var');
    FWHM = 0; % no smoothing
end
hemi = {'lh' 'rh'};
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
FS_atlas = fullfile(SUBJECTS_DIR,'cvs_avg35_inMNI152/mri/T1.nii.gz');
FSL_DIR = getenv('FSL_DIR');
FSL_atlas = fullfile(FSL_DIR,'data/standard/MNI152_T1_2mm.nii.gz');
%% Create registration between FSL 2mm MNI and Freesurfer 1mm MNI space
disp('Registering atlases...');
[~,~] = system(['tkregister2 --mov ' FSL_atlas ' --targ ' FS_atlas ' --regheader ' ...
    ' --reg ~/tmp.dat --noedit']);
%% Project from FSL 2mm MNI space to Freesurfer 1mm MNI space
[~,~] = system(['mri_vol2vol --mov ' input_vol ' --targ ' ...
    fullfile(SUBJECTS_DIR,'cvs_avg35_inMNI152/mri/T1.nii.gz') ' --o ' ...
    '~/tmp.nii.gz --reg ~/tmp.dat']);
%% Project to Freesurfer cvs_avg35_inMNI152 surface
for hh = 1:length(hemi)
    disp(['Projecting to ' hemi{hh} ' surface...']);
    out_surf = fullfile(out_path,[hemi{hh} '.' outname '.nii.gz']);
 [~,~] = system(['mri_vol2surf --src ~/tmp.nii.gz --regheader ' ...
     'cvs_avg35_inMNI152 --hemi ' hemi{hh} ' --surf-fwhm ' num2str(FWHM) ...
     ' --out ' out_surf ' --projfrac 0.5']);
end
%% Clean up files
!rm ~/tmp.dat ~/tmp.nii.gz*
disp('done.');
