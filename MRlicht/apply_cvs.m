function apply_cvs(in_vol,out_vol,subject_name,SUBJECTS_DIR)

% Apply cvs registration created using mri_cvs_register
%
%   Assumes the subject was registered using:
%
%   mri_cvs_register --mov <subject_name> --template cvs_avg35_inMNI152
%
%   Usage:
%
%   apply_cvs(in_vol,out_vol,subject_name,SUBJECTS_DIR)
%
%   Written by Andrew S Bock May 2015

%   07/11/15 kmo, asb    Removed some stuff that were buggy in Maria's
%                        pipeline (and apparently unnecessary).

%% Set Defaults
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR=getenv('SUBJECTS_DIR');
end
cvsmorphfile = fullfile(SUBJECTS_DIR,subject_name,...
    'cvs/final_CVSmorph_tocvs_avg35_inMNI152.m3z');
targfile = fullfile(SUBJECTS_DIR,'cvs_avg35_inMNI152/mri/T1.mgz');
%% Apply cvs registration
disp(['Applying cvs registration to ' in_vol '...']);
system(['mri_vol2vol --targ ' targfile ...
    ' --noDefM3zPath --m3z ' cvsmorphfile ...
    ' --no-save-reg --mov ' in_vol ' --o ' out_vol ...
    ' --interp nearest']);
disp('done.');