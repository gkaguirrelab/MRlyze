function calc_M_vs_P_voxels(session_dir,func,gfeatNum)

%   This will create a volume with values than range between 0 (P) and 1
%   (M)
%
%   Usage:
%   calc_M_vs_P_voxels(session_dir,func,gfeatNum)
%
%   Written by Andrew S Bock Nov 2015

%% Find bold directories
b = find_bold(session_dir);
featDir = fullfile(session_dir,b{gfeatNum},[func '.gfeat']);
%% Get voxel 'cope_1_2_prop.nii.gz' values
cope1 = load_nifti(fullfile(featDir,'cope1.feat','stats','cope1.psc.nii.gz'));
cope2 = load_nifti(fullfile(featDir,'cope2.feat','stats','cope1.psc.nii.gz'));
cope3 = load_nifti(fullfile(featDir,'cope3.feat','stats','cope1.psc.nii.gz'));
cope4 = load_nifti(fullfile(featDir,'cope4.feat','stats','cope1.psc.nii.gz'));
copesum = cope1;
% M vs P
copesum.vol = abs(cope1.vol)./(abs(cope1.vol)+abs(cope2.vol));
save_nifti(copesum,fullfile(session_dir,[func '.M_vs_P.nii.gz']));
% M>P_vs_P<M
copesum.vol = abs(cope3.vol)./(abs(cope3.vol)+abs(cope4.vol));
save_nifti(copesum,fullfile(session_dir,[func '.MP_vs_PM.nii.gz']));
