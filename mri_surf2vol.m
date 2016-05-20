function mri_surf2vol(subject_name,invol,outvol,hemi,trgvol)

% Projects from the freesurfer cortical surface to volume
%
%   Usage:
%   mri_surf2vol(subject_name,invol,outvol,hemi,trgvol,SUBJECTS_DIR)
%
%   Written by Andrew S Bock Nov 2015

%% Set defaults
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
if ~exist('trgvol','var') || isempty(trgvol)
    trgvol = fullfile(SUBJECTS_DIR,subject_name,'mri','orig.mgz');
end
%% Project from subject surface space to volume
system(['mri_surf2vol --surfval ' invol ' --hemi ' hemi ' --fillribbon --identity ' ...
    subject_name ' --template ' trgvol ' --o ' outvol]);
