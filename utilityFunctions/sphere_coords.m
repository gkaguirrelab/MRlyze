function [SphCoords] = sphere_coords(mask,radius,voxsize)

%   Takes in a binary 3D mask volume, calculates the surrounding voxels
%   within a given radius for each mask voxel (i.e. = 1) in the volume.
%
%   Note: the indices that result are in the volume, but they are organized
%   by the mask. I.e. the size of the output structure is based on the size
%   of the mask, but the sphere coordinates for each mask voxel in the
%   output structure are in terms of the entire volume.
%
%   Usage:
%   [out] = vol_sphere(mask,radius,voxsize)
%
%   Defaults:
%   mask - no default, must be a binary volume
%   radius - 20 (mm)
%   voxsize - 1.5 (mm)
%
%   Output:
%   out - structure with values for each mask voxel (v), with the following
%   fields for each voxel:
%       mask_coords - coordinates of voxels in the mask (Vx x 3)
%       vol_idx - volume index of voxels in the sphere
%       neighb_counts - number of voxels in the sphere
%
%   Modeled after the 'adj_sphere' function in the Princeton MVPA toolbox
%   and the Anaticor toolbox.
%
%   Note: I choose 20 mm as the default, as I found subcortical structures (e.g. SC,
%   IC) require a slightly larger sphere than the typical 15mm in
%   Anaticor.
%
%   Written by Andrew S Bock Sept 2014

%% Set up initial variables
if ~exist('mask','var')
    error('No mask volume')
end
if ~exist('radius','var')
    radius = 20;
end
if ~exist('voxsize','var')
    voxsize = 1.5; % assumes isotropic voxel size
end
%% Create sphere
% offsets is the number of voxels, so we scale by voxel size
disp('Creating initial sphere');
ceil_radius = ceil(radius/voxsize); % allows for non-interger radii
offsets = [];
for x = -1*ceil_radius:ceil_radius
    for y = -1*ceil_radius:ceil_radius
        for z = -1*ceil_radius:ceil_radius
            if sqrt(x^2+y^2+z^2) <= radius
                % append the current offset-triplet to the next row
                offsets(end+1,:) = [x y z];
            end
        end
    end
end
nvox_in_sphere = size(offsets,1);
%% Find the x,y,z coordinates of the mask, in terms of the overall volume.
mask_dims = size(mask);
[mask_coords(:,1),mask_coords(:,2),mask_coords(:,3)] = ind2sub(mask_dims,find(mask));
SphCoords.mask_coords = mask_coords;
nvox_in_mask = size(mask_coords,1);
% Initialize various variables
disp('Initializing Sphere Coordinate variables');
SphCoords.vol_idx = zeros(nvox_in_mask, nvox_in_sphere, 'single');
SphCoords.neighb_counts = zeros(nvox_in_mask,1,'single');
%% Place sphere around each voxel
progBar = ProgressBar(nvox_in_mask, ...
    'Finding surrounding voxels in sphere for each voxel in mask');
for v = 1:nvox_in_mask
    vox_coords_tiled = repmat(mask_coords(v,:), nvox_in_sphere, 1);
    sphere_coords = vox_coords_tiled + offsets;
    % Remove voxels outside of the volume
    badvoxels = sphere_coords(:,1)<1 | ...
        sphere_coords(:,2)<1 | ...
        sphere_coords(:,3)<1 | ...
        sphere_coords(:,1)>mask_dims(1) | ...
        sphere_coords(:,2)>mask_dims(2) | ...
        sphere_coords(:,3)>mask_dims(3);
    sphere_coords(badvoxels,:) = [];
    % Get the idicies of voxels within the sphere
    sphere_idx_vol = sub2ind(mask_dims, sphere_coords(:,1), ...
        sphere_coords(:,2), sphere_coords(:,3));
    % Number of voxels in sphere
    nvox_remaining_in_sphere = length(sphere_idx_vol);
    % Update matrix with index values in the overall volume
    SphCoords.vol_idx(v, 1:nvox_remaining_in_sphere) = sphere_idx_vol;
    % Update vector saving the number of voxels in the sphere
    SphCoords.neighb_counts(v) = nvox_remaining_in_sphere;
    if ~mod(v,round(nvox_in_mask/20));progBar(v);end
end