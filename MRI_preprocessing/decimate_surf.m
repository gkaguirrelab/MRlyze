function decimate_surf(subject_name,hemi,src_surf,trg_surf,in_vol,out_vol)

% Decimates surface image (3D or 4D).  
%
%   Usage:
%   decimate_surf(subject_name,hemi,src_surf,trg_surf,in_vol,out_vol)
%
%   Example:
%   subject_name = 'ASB';
%   hemi = 'lh';
%   src_surf = 'inflated';
%   trg_surf = '0.1.inflated';
%   in_vol = '/path/to/in/surf/file.nii.gz';
%   out_vol = '/path/to/out/surf/file.nii.gz';
%   decimate_surf(subject_name,hemi,src_surf,trg_surf,in_vol,out_vol);
%
%   Assumes that the following has been done in terminal (must be Linux!)
%
%   cd $SUBJECTS_DIR/<subject_name>/surf
%   mris_decimate -d 0.1 ./<hemi>.<surf> ./<hemi>.0.1.<surf>
%
%   e.g.
%   mris_decimate -d 0.1 ./lh.inflated ./lh.0.1.inflated
%   mris_decimate -d 0.1 ./rh.inflated ./rh.0.1.inflated
%
%   Written by Andrew S Bock Oct 2015

%% Decimate input
disp('Decimating...');
[sind] = find_closest_verts(subject_name,hemi,src_surf,trg_surf);
project_2_decimate(in_vol,out_vol,sind);
disp('done.');
