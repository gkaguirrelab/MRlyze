function apply_cvs_inverse(in_vol,out_vol,ref_vol,subject_name,SUBJECTS_DIR)

% Apply cvs registration created using mri_cvs_register, in the inverse
% direction
%
%   Usage:
%
%   apply_cvs_inverse(in_vol,out_vol,ref_vol,subject_name,SUBJECTS_DIR)
%
%   Assumes the subject was registered using:
%
%   mri_cvs_register --mov <subject_name> --template cvs_avg35_inMNI152
%
%   Written by Andrew S Bock Jul 2015

%% Set Defaults
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR=getenv('SUBJECTS_DIR');
end
cvsmorphfile = fullfile(SUBJECTS_DIR,subject_name,...
    'cvs/final_CVSmorph_tocvs_avg35_inMNI152.m3z');
%% Apply cvs registration
disp(['Applying cvs registration to ' in_vol '...']);
system(['mri_vol2vol --targ ' in_vol ...
    ' --noDefM3zPath --m3z ' cvsmorphfile ...
    ' --no-save-reg --mov ' ref_vol ' --o ' out_vol ...
    ' --inv-morph --interp nearest']);
disp('done.');