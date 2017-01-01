function project_CF_cvsMNI(session_dir,subject_name,map_type,maps,thresh)

% Averages CF volume maps in cvsMNI space (freesurfer) using the
% 'average_in_cvsMNI' function.
%
%   Usage: 
%   project_CF_cvsMNI(session_dir,subject_name,maps,thresh)
%
%   Defaults:
%   maps = {'var' 'varsig1' 'varsig2' 'varsig3' 'varsig4' 'varecc' 'varpol'};
%   thresh = 0.125; % variance explained
%
%   Written by Andrew S Bock May 2015
%% Set defaults
if ~exist('maps','var')
    maps = {'var' 'varsig1' 'varsig2' 'varsig3' 'varsig4' 'varecc' 'varpol'};
end
if ~exist('thresh','var')
    thresh = 0.125; % variance explained
end
%% Average in cvsMNI space
for m = 1:length(maps)
    lh_vol_in = fullfile(session_dir,['lh.volume.prf_V1.' map_type '.avg.' maps{m} '.nii.gz']);
    rh_vol_in = fullfile(session_dir,['rh.volume.prf_V1.' map_type '.avg.' maps{m} '.nii.gz']);
    mh_vol_out = fullfile(session_dir,['mh.volume.prf_V1.' map_type '.avg.' maps{m} '.mni.nii.gz']);
    average_in_cvsMNI(session_dir,subject_name,maps{m},lh_vol_in,rh_vol_in,mh_vol_out)
end
%% Threshold resulting maps
varvol = load_nifti(fullfile(session_dir,['mh.volume.prf_V1.' map_type '.avg.var.mni.nii.gz']));
varthresh = varvol.vol>=thresh;
disp(['Thresholding maps at ' num2str(thresh) ' variance explained...']);
for m = 1:length(maps)
    mh_vol_in = fullfile(session_dir,['mh.volume.prf_V1.' map_type '.avg.' maps{m} '.mni.nii.gz']);
    mh_vol_out = fullfile(session_dir,['mh.volume.prf_V1.' map_type '.avg.' maps{m} '.mni.thresh.nii.gz']);
    nii = load_nifti(mh_vol_in);
    nii.vol(~varthresh) = nan;
    save_nifti(nii,mh_vol_out);
end
disp('done.');