function decimate_bold(session_dir,subject_name,func,src_surf,trg_surf)

% Decimates bold 4D volumes.  Assumes that the following has been done in
% terminal (must be Linux!).
%
%   cd $SUBJECTS_DIR/<subject_name>/surf
%   mris_decimate -d 0.1 ./lh.inflated ./lh.0.1.inflated
%   mris_decimate -d 0.1 ./rh.inflated ./rh.0.1.inflated
%
%   Usage:
%   decimate_bold(session_dir,subject_name)
%
%   Written by Andrew S Bock Sep 2015

%% Set defaults
bolddirs = find_bold(session_dir);
if ~exist('func','var')
   func = 's5.dbrf.tf'; 
end
if ~exist('src_surf','var')
    src_surf = 'inflated';
end
if ~exist('trg_surf','var')
    trg_surf = '0.1.inflated';
end
hemis = {'lh' 'rh'};
%% Decimate bold timecourses
disp('Decimating bold timecourses...');
for i = 1:length(bolddirs)
    for hh = 1:2
        hemi = hemis{hh};
        [sind] = find_closest_verts(subject_name,hemi,src_surf,trg_surf);
        in_name = fullfile(session_dir,bolddirs{i},[func '.surf.' hemi '.nii.gz']);
        out_name = fullfile(session_dir,bolddirs{i},[func '.surf.' trg_surf '.' hemi '.nii.gz']);
        project_2_decimate(in_name,out_name,sind);
    end
end