function mri_surf2vol(subject_name,invol,outvol,hemi,templatevol,SUBJECTS_DIR)

% Projects the pRF maps ('ecc','pol','areas') from surface to volume
%
%   Usage:
%   mri_surf2vol(session_dir,subject_name,invol,outvol,hemi,templatevol,SUBJECTS_DIR)
%
%   Written by Andrew S Bock Nov 2015

%% Set defaults
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
if ~exist('templatevol','var')
    templatevol = fullfile(SUBJECTS_DIR,subject_name,'mri','T1.mgz');
end
%% Project from subject surface space to volume
system(['mri_surf2vol --surfval ' invol ' --hemi ' hemi ' --fillribbon --identity ' ...
    subject_name ' --template ' templatevol ' --o ' outvol]);
