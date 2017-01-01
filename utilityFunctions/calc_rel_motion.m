function calc_rel_motion(rawvol,mcDir)

% Calculates the relative motion based on the Euclidean distance of each
% voxels position before and after motion correction, for each TR.
%
%   Usage:
%   calc_rel_motion(rawvol,mcDir)
%
%   Inputs:
%   rawvol = raw volume, prior to any motion correctoin
%   mcDir = motion correction directory (e.g.
%   <featDir>/mc/prefiltered_func_data_mcf.mat)
%
%   This follows:
%
%   Satterthwaite, T. D., Elliott, M. A., Gerraty, R. T., Ruparel, K.,
%   Loughead, J., Calkins, M. E., ... & Wolf, D. H. (2013). An improved
%   framework for confound regression and filtering for control of motion
%   artifact in the preprocessing of resting-state functional connectivity
%   data. Neuroimage, 64, 240-256.
%
%   From http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html
%
%   The voxel index coordinates (i.e., the array indexes) are referred to
%   as (i,j,k), with valid ranges:
%       i = 0 .. dim[1]-1
%       j = 0 .. dim[2]-1  (if dim[0] >= 2)
%       k = 0 .. dim[3]-1  (if dim[0] >= 3)
%
%       x = srow_x[0] * i + srow_x[1] * j + srow_x[2] * k + srow_x[3]
%       y = srow_y[0] * i + srow_y[1] * j + srow_y[2] * k + srow_y[3]
%       z = srow_z[0] * i + srow_z[1] * j + srow_z[2] * k + srow_z[3]
%
%   Written by Andrew S Bock Sep 2015
%% Load unprocessed functional volume, pull out dims and sform
raw = load_nifti(rawvol);
dims = size(raw.vol);
sform = raw.sform;
%% Load in motion correction matricies for each TR
m = listdir(fullfile(mcDir,'MAT_*'),'files');
for i = 1:length(m)
    motMat(i,:,:) = load(m{i});
end
%% Get voxel scanner coordinates
for x = 1:dims(1)
    for y = 1:dims(2)
        for z = 1:dims(3)
            % set i,j,k location for voxel
            vi = x-1; vj = y-1; vk = z-1;
            vxCoords(x,y,z,1) = sform(1,1)*vi + sform(1,2)*vj + sform(1,3)*vk + sform(1,4);
            vxCoords(x,y,z,2) = sform(2,1)*vi + sform(2,2)*vj + sform(2,3)*vk + sform(2,4);
            vxCoords(x,y,z,3) = sform(3,1)*vi + sform(3,2)*vj + sform(3,3)*vk + sform(3,4);
            vxCoords(x,y,z,4) = 1;
        end
    end
end
%% Apply affine transform
relMot = nan(dims);
for i = 1%:dims(4)
    A = squeeze(motMat(i,:,:));
    for x = 50%dims(1)
        for y = 1:dims(2)
            for z = 1:dims(3)
                tmpVxCoords = (squeeze(vxCoords(x,y,z,:)))';
                newvxCorrds = tmpVxCoords*A';
                relMot(x,y,z,i) = sqrt( ...
                    (tmpVxCoords(1) - newvxCorrds(1))^2 + ...
                    (tmpVxCoords(2) - newvxCorrds(2))^2 + ...
                    (tmpVxCoords(3) - newvxCorrds(3))^2 ...
                    );
            end
        end
    end
end
%%
% To express instantaneous head motion as a scalar quantity we used the empirical formula, FDi = |?dix| + |?diy| + |?diz| + |??i| + |??i| + |??i|, where ?dix = d(i ? 1)x ? dix, and similarly for the other rigid body parameters [dix diy diz ?i ?i ?i]. 
% Rotational displacements were converted from degrees to millimeters by calculating displacement on the surface of a sphere of radius 50 mm, which is approximately the mean distance from the cerebral cortex to the center of the head.