function make_SC_sphere(session_dir,subject_name,center_voxs,func,SUBJECTS_DIR)

% Creates spheres around center voxels for the left and right superior
% colliculi, for use as masks in later analyses.
%
%   Usage:
%   make_SC_sphere(session_dir,subject_name,center_voxs,func,SUBJECTS_DIR)
%
%   Default numVox values of (500mm3) are initally based on:
%
%       SC (~86mm3)
%       Petit, L., & Beauchamp, M. S. (2003). Journal of neurophysiology, 89(5), 2516-2527.
%       Schneider, K. A., & Kastner, S. (2005). Journal of Neurophysiology, 94(4), 2491-2503.
%
%   Note: the volume estimate of SC (500mm3) is ~5x the reported mean SC volume,
%   thus it is an over-estimate.
%
%   Written by Andrew S Bock Sep 2015

%% set defaults
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
if ~exist('func','var')
    func = 'brf'; % functional data file used for registration
end
voxsize = 1; % assumes 1mm isotropic voxel size
input_vol = fullfile(session_dir,['mean_' func '.anat.nii.gz']);
output_vol{1} = fullfile(session_dir,'lh.SC.nii.gz');
output_vol{2} = fullfile(session_dir,'rh.SC.nii.gz');
ROIsize = 500; % default for SC, mm^3
radius = ceil( ((ROIsize)/(4/3*pi)) ^ (1/3) ); % find radius, given ROIsize volume
%% Find bold run directories
d = find_bold(session_dir);
%% Project functional runs to anatomical space
if ~exist(input_vol,'file');
    for i = 1:length(d);
        runNum = i;
        project_func2anat(session_dir,subject_name,runNum,func,SUBJECTS_DIR)
    end
    % Create average functional image in anatomical space
    tmp = [];
    for i = 1:length(d);
        runNum = i;
        tmpvol = load_nifti(fullfile(session_dir,d{runNum},[func '.anat.nii.gz']));
        tmp(i,:,:,:) = tmpvol.vol;
    end
    tmpvol.vol = squeeze(mean(tmp));
    save_nifti(tmpvol,input_vol);
end
%% Load the data
tmp = load_nifti(input_vol);
%% Create sphere
ceil_radius = ceil(radius/voxsize); % allows for non-interger radii
offsets = [];
for x = -1*ceil_radius:ceil_radius
    for y = -1*ceil_radius:ceil_radius
        for z = -1*ceil_radius:ceil_radius
            if sqrt(x^2+y^2+z^2) <= ceil_radius
                % append the current offset-triplet to the next row
                offsets(end+1,:) = [x y z];
            end
        end
    end
end
nvox_in_sphere = size(offsets,1);
%% Create sphere around SC in cvs space
for i = 1:length(center_voxs)
    center_vox = center_voxs{i};
    vox_coords_tiled = repmat(center_vox, nvox_in_sphere, 1);
    % create sphere around GM voxel
    tmp_coords = vox_coords_tiled + offsets;
    % Remove voxels outside of the volume
    mask_dims = size(tmp.vol);
    badvoxels = tmp_coords(:,1)<1 | ...
        tmp_coords(:,2)<1 | ...
        tmp_coords(:,3)<1 | ...
        tmp_coords(:,1)>mask_dims(1) | ...
        tmp_coords(:,2)>mask_dims(2) | ...
        tmp_coords(:,3)>mask_dims(3);
    sphere_coords = tmp_coords(~badvoxels,:);
    tmp.vol = zeros(size(tmp.vol));
    for j = 1:size(sphere_coords,1);
        tmp.vol(sphere_coords(j,1),sphere_coords(j,2),sphere_coords(j,3)) = 1;
    end
    save_nifti(tmp,output_vol{i});
end
%% Combine across hemispheres
lh = load_nifti(output_vol{1});
rh = load_nifti(output_vol{2});
mh = lh;
mh.vol = lh.vol + rh.vol;
mh.vol(mh.vol>0) = 1;
save_nifti(mh,fullfile(session_dir,'mh.SC.nii.gz'));
%% Deprecated
% %% Project to subject space
% for i = 1:length(center_voxs)
%     in_vol = tmp_output{i};
%     out_vol = output_vol{i};
%     ref_vol = fullfile(SUBJECTS_DIR,subject_name,'mri','T1.mgz');
%     apply_cvs_inverse(in_vol,out_vol,ref_vol,subject_name,SUBJECTS_DIR)
% end
% %% Remove temporary volumes
% for i = 1:length(center_voxs)
%     [~,~] = system(['rm ' tmp_output{i}]);
% end