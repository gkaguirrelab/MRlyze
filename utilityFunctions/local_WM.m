function [wm_tc] = local_WM(brain,WM,tc,voxsize,radius)

% Takes in brain and white matter binary volumes, as well as their
%   corresponding timecourses. Modeled after the AFNI Anaticor pipeline.
%
% Usage:
%   [wm_tc] = local_WM(brain,WM,tc,voxsize,radius)
%
% Defaults:
%   radius - 15 (mm)
%
% Output:
%   The average local WM timecourse for each brain voxel
%
%   Note: Sphere radius default = 15mm based on Anaticor.
%
%   Note: this function uses parpool. If you have trouble starting a
%   parallel pool, due to java exception errors, try adding
%
%   127.0.0.1 <hostname>
%
%   to /etc/hosts.  e.g. 127.0.0.1 mac-pro.uphs.upenn.edu
%
%   Written by Andrew S Bock Sept 2014

%% Set up initial variables
if ~exist('radius','var')
    radius = 15;
end
%% Create sphere
% offsets is the number of voxels, so we scale by voxel size
disp(['Creating initial ' num2str(radius) 'mm sphere...']);
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
disp('done.');
%% Find the x,y,z coordinates of the GM mask, WM indices, and flat timecourses
% Brain coords and mask
mask_dims = size(brain);
brainind = find(brain);
[brain_coords(:,1),brain_coords(:,2),brain_coords(:,3)] = ind2sub(mask_dims,brainind);
nvox_in_mask = size(brain_coords,1);
% WM indices
WMind = logical(reshape(WM,[1,mask_dims(1)*mask_dims(2)*mask_dims(3)]));
% Find voxels with zero variance (i.e. flat timecourses)
tcvar = var(tc);
flatind = tcvar == 0; 
% combine non-WM and flat timecourses
badind = ~WMind | flatind; 
%% Place sphere around each voxel
%progBar = ProgressBar(nvox_in_mask,'Finding surrounding voxels in sphere for each voxel in mask');
tmp_tc = zeros(size(tc,1),nvox_in_mask);
wm_tc = zeros(size(tc));
% Slice for parfor
mask_dims_one = mask_dims(1);
mask_dims_two = mask_dims(2);
mask_dims_three = mask_dims(3);
poolobj = gcp; % Gets current pool, and if no pool, creates a new one
disp('Calculating local white matter timecourse for each voxel...');
tstart = clock;
ProgressBar_parfor(tstart);
%progBar=ProgressBar(nvox_in_mask,'foo');
parfor v = 1:nvox_in_mask
    % tile GM voxel X,Y,Z
    vox_coords_tiled = repmat(brain_coords(v,:), nvox_in_sphere, 1);
    % create sphere around GM voxel
    tmp_coords = vox_coords_tiled + offsets;
    % Remove voxels outside of the volume
    badvoxels = tmp_coords(:,1)<1 | ...
        tmp_coords(:,2)<1 | ...
        tmp_coords(:,3)<1 | ...
        tmp_coords(:,1)>mask_dims_one | ...
        tmp_coords(:,2)>mask_dims_two | ...
        tmp_coords(:,3)>mask_dims_three;
    sphere_coords = tmp_coords(~badvoxels,:);
    % Convert sphere_coords from X,Y,Z to ind
    sphere_coords_ind = sub2ind(mask_dims,sphere_coords(:,1),sphere_coords(:,2),sphere_coords(:,3));
    % Initialize logical for sphere
    sphere_ind = false(mask_dims(1)*mask_dims(2)*mask_dims(3),1);
    % Set sphere coordinates to true
    sphere_ind(sphere_coords_ind) = 1;
    sphere_ind(badind) = 0; % Set non-WM and flat timecourses to zero
    % Get the average timecourse of the WM voxels in the sphere
    tmp_tc(:,v) = mean(tc(:,sphere_ind),2);
    %     %if ~mod(v,round(nvox_in_mask/1000));progBar(v);end
    %progBar(v);
    ProgressBar_parfor(tstart,v,nvox_in_mask);
end
wm_tc(:,brainind) = tmp_tc;
wm_tc(isnan(wm_tc)) = 0;
disp('done.');
ProgressBar_parfor(tstart,'clean'); % Clean up files and display total loop time
delete(poolobj); % close parpool