function [arclen] = spherical_distance(subject_name,vtx,verts,hemi)

% Finds the arc length between a point and all other points on the cortical
% surface, using the subject's 'sphere' surface.
%
%   Usage:
%   [arclen] = spherical_distance(subject_name,vtx,verts,hemi)
%
%   assumes 'verts' was created as follows:
%       anatdatadir = fullfile(SUBJECTS_DIR,subject_name);
%       [verts] = freesurfer_read_surf(fullfile(anatdatadir,'surf',[hemi '.sphere']));
%
%   Written by Andrew S Bock Jul 2015

%% Set defaults
if ~exist('hemi','var')
    hemi = 'lh';
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
anatdatadir = fullfile(SUBJECTS_DIR,subject_name);
%% Calculate arc length
[theta_sph,phi_sph,rho_sph] = cart2sph(verts(:,1),verts(:,2),verts(:,3));
tmparc = rho_sph(vtx)*distance(theta_sph(vtx),phi_sph(vtx),theta_sph(:),phi_sph(:));
%% Scale arc length by surface area
% We do this because the sphere surface assumes a radius of 100mm, but
%   individual brains vary in size.
% Load stat file
statfile = fullfile(anatdatadir,'stats',[hemi '.aparc.stats']);
fid   = fopen(statfile);
A     = fread(fid,'char');
Achar = char(A');
fclose(fid);
% Find surface area
idx1    = strfind(Achar, 'White Surface Total Area');
idx2    = strfind(Achar, ', mm^2');
surface_area = str2double(Achar(idx1+25:idx2-1));
sphere_area = (4*pi*rho_sph(vtx)^2);
% Scale by the sqrt of the surface area ratio
scale_factor = sqrt(surface_area / sphere_area);
arclen = tmparc * scale_factor;