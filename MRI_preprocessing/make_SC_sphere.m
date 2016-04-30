function make_SC_sphere(inVol,outVols,center_voxs,SCvolume)

% Creates spheres around center voxels for the left and right superior
% colliculi, for use as masks in later analyses.
%
%   Usage:
%   make_SC_sphere(inVol,outVols,center_voxs,SCvolume)
%
%   Inputs:
%   inVol           = input anatomical volume (1x1x1mm)
%   center_voxs{1}  = center voxels (x,y,z) for lh SC
%   center_voxs{2}  = center voxels (x,y,z) for rh SC
%   center_voxs{3}  = center voxels (x,y,z) for mh SC
%
%   Default:
%   SCvolume        = 500; 
%
%       SC (~86mm3)
%       Petit, L., & Beauchamp, M. S. (2003). Journal of neurophysiology, 89(5), 2516-2527.
%       Schneider, K. A., & Kastner, S. (2005). Journal of Neurophysiology, 94(4), 2491-2503.
%
%   Note: the volume estimate of SC (500mm3) is ~6x the reported mean SC volume,
%   thus it is an over-estimate.
%
%   Written by Andrew S Bock Sep 2015

%% set defaults
if ~exist('SCvolume','var')
    SCvolume = 500; % default for SC, mm^3
end
radius = ceil( ((SCvolume)/(4/3*pi)) ^ (1/3) ); % find radius, given SCvolume volume
%% Load input volume
anat = load_nifti(inVol);
%% Create sphere
ceil_radius = ceil(radius); % allows for non-interger radii
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
    mask_dims = size(anat.vol);
    badvoxels = tmp_coords(:,1)<1 | ...
        tmp_coords(:,2)<1 | ...
        tmp_coords(:,3)<1 | ...
        tmp_coords(:,1)>mask_dims(1) | ...
        tmp_coords(:,2)>mask_dims(2) | ...
        tmp_coords(:,3)>mask_dims(3);
    sphere_coords = tmp_coords(~badvoxels,:);
    anat.vol = zeros(size(anat.vol));
    for j = 1:size(sphere_coords,1);
        anat.vol(sphere_coords(j,1),sphere_coords(j,2),sphere_coords(j,3)) = 1;
    end
    save_nifti(anat,outVols{i});
end
%% Combine across hemispheres
if length(outVols) == 3
    lh = load_nifti(outVols{1});
    rh = load_nifti(outVols{2});
    mh = lh;
    mh.vol = lh.vol + rh.vol;
    mh.vol(mh.vol>0) = 1;
end
save_nifti(mh,outVols{3});