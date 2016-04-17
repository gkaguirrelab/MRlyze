function create_subcortical_CF_masks(session_dir,feat_session_dir,ROI,hemi,center_vox)

% Creates subcortical maps using the correlation values resulting from the
% CF pipeline
%
%   Usage:
%   create_subcortical_CF_maps(session_dir,feat_session_dir,ROI,hemi,center_vox,srcROI,trgROI,map_type)
%   ROI = 'LGN';
%   center_vox = [149 154 121];
%
%   Written by Andrew S Bock Jun 2015

%% Find the gfeat dir
d = listdir(fullfile(feat_session_dir,'*BOLD_*'),'dirs');
if isempty(d)
    d = listdir(fullfile(feat_session_dir,'*bold_*'),'dirs');
end
if isempty(d)
    d = listdir(fullfile(feat_session_dir,'*EPI_*'),'dirs');
end
if isempty(d)
    d = listdir(fullfile(feat_session_dir,'*ep2d*'),'dirs');
end
if isempty(d)
    d = listdir(fullfile(feat_session_dir,'RUN*'),'dirs');
end
%% Create mask
if strcmp(hemi,'rh')
    stat_vol = fullfile(feat_session_dir,d{1},'sdbrf.tf.gfeat/cope3.feat/stats/zstat1.nii.gz');
else
    stat_vol = fullfile(feat_session_dir,d{1},'sdbrf.tf.gfeat/cope4.feat/stats/zstat1.nii.gz');
end
mask_vol = fullfile(session_dir,[hemi '.' ROI '.mask.nii.gz']);
find_ROI_voxels(stat_vol,mask_vol,center_vox,ROI);
mask = load_nifti(mask_vol);
stats = load_nifti(stat_vol);
stats.vol = stats.vol .* mask.vol;
save_nifti(stats,fullfile(session_dir,[hemi '.' ROI '.mask.zstat1.nii.gz']));
