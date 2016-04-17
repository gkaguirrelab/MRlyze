function [sind] = find_closest_verts(subject_name,hemi,src_surf,trg_surf,SUBJECTS_DIR)

% Finds the closest vertices, useful for a decimated surface. Writes out a
% decimated *.curv file as well.
%
% Usage:
%   [sind] = find_closest_verts(subject_name,hemi,src_surf,trg_surf,SUBJECTS_DIR)
%
%   output:
%   sind - closest vertices in the source surface for the decimate surface
%
%   example:
%   find_closest_verts('A102714B','lh','sphere','0.1.sphere');
%
%   Written by Andrew S Bock Aug 2015


%% set defaults
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
if ~exist('src_surf','var')
    src_surf = 'sphere';
end
if ~exist('trg_surf','var')
    trg_surf = '0.1.sphere';
end
%% Load in the surfaces
[svert,~] = freesurfer_read_surf(...
    fullfile(SUBJECTS_DIR,subject_name,'surf',[hemi '.' src_surf]));
[tvert,tface] = freesurfer_read_surf(...
    fullfile(SUBJECTS_DIR,subject_name,'surf',[hemi '.' trg_surf]));
%% Find closest vertex in target surface
sind = nan(length(tvert),1);
disp('calculating distance between vertices...');
for i = 1:length(tvert)
    d = sqrt(...
        (tvert(i,1) - svert(:,1)).^2 + ...
        (tvert(i,2) - svert(:,2)).^2 + ...
        (tvert(i,3) - svert(:,3)).^2);
    [~,sind(i)] = min(d);
end
disp('done.');
%% Write out curvature file
[icurv,~] = freesurfer_read_curv(...
    fullfile(SUBJECTS_DIR,subject_name,'surf',[hemi '.curv']));
out_curv = fullfile(SUBJECTS_DIR,subject_name,'surf',[hemi '.' trg_surf '.curv']);
ocurv = icurv(sind);
if ~exist(out_curv,'file');
    write_curv(out_curv, ocurv, length(tface));
else
    disp(['WARNING: ' out_curv ' already exists!']);
end
