function [ROI_vox,ROI_vals] = find_ROI_voxels(input_vol,output_vol,center_vox,ROI,voxsize)

% Finds a specified number of voxels within a sphere surrounding a center
% voxel. Voxels are chosen in a decending order (i.e. find the highest
% values).
%
%   Usage:
%   [ROI_vox,ROI_vals] = find_ROI_voxels(input_vol,output_vol,center_vox,ROI,voxsize)
%
%   Default numVox values of LGN (200mm3) and SC (75mm3) are based on:
%
%       LGN (~260mm3)
%       Korsholm, K. et ak (2007).Brain, 130(5), 1244-1253.
%       O'Connor, D. H., ... Kastner, S. (2002). Nature neuroscience, 5(11), 1203-1209.
%       Chen, W., ... Ugurbil, K. (1999). Proceedings of the National Academy of Sciences, 96(5), 2430-2434.
%
%       SC (~86mm3)
%       Petit, L., & Beauchamp, M. S. (2003). Journal of neurophysiology, 89(5), 2516-2527.
%       Schneider, K. A., & Kastner, S. (2005). Journal of Neurophysiology, 94(4), 2491-2503.
%
%
%   Note: the number of voxels is slightly smaller (-1 SEM) than the reported mean
%   LGN and SC volumes, thus it is a conservative estimate.
%
%   Written by Andrew S Bock Jun 2015

%% Set defaults
if strcmp(ROI,'LGN')
    ROIsize = 200; % default for LGN
elseif strcmp(ROI,'SC')
    ROIsize = 75; % default for SC
end
if ~exist('voxsize','var')
    voxsize = 1.5; % assumes isotropic voxel size
end
numVox = ceil(ROIsize/voxsize);
radius = ceil( ((3*ROIsize)/(4/3*pi)) ^ (1/3) ); % triple the size of the expected ROIsize
%% Load the data
disp('Finding ROI voxels...');
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
%% tile GM voxel X,Y,Z
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
%% Find stat values in sphere
stat_vals = zeros(size(sphere_coords,1),1);
for i = 1:size(sphere_coords,1);
    stat_vals(i) = tmp.vol(sphere_coords(i,1),sphere_coords(i,2),sphere_coords(i,3));
end
%% Find max numVox
[Y,I] = sort(stat_vals);
ROI_vox = sphere_coords(I(end-numVox+1:end),:);
ROI_vals = Y(end-numVox+1:end);
%% Save mask
tmp.vol = zeros(size(tmp.vol));
for i = 1:numVox
    tmp.vol(ROI_vox(i,1),ROI_vox(i,2),ROI_vox(i,3)) = 1;
end
save_nifti(tmp,output_vol);
disp('done.');
