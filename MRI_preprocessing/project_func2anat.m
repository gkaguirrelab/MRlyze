function project_func2anat(session_dir,subject_name,runNum,func,SUBJECTS_DIR)
% Projects a functional file (single TR) to anatomical space
%
%   Usage:
%   project_func2anat(session_dir,subject_name,runNum,func,SUBJECTS_DIR)
%
%   Written by Andrew S Bock Oct 2015

%% Set defaults
if ~exist('runNum','var')
    runNum = 1; % functional data file used for registration
end
if ~exist('func','var')
    func = 'brf'; % functional data file used for registration
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
%% Find bold run directories
d = find_bold(session_dir);
%% Pull out first TR
invol = fullfile(session_dir,d{runNum},[func '.nii.gz']);
outvol = fullfile(session_dir,d{runNum},'tmp.nii.gz');
[~,~] = system(['fslroi ' invol ' ' outvol ' 0 1']); % create a 3D volume
%% Project functional volume to anatomical space
invol = fullfile(session_dir,d{runNum},'tmp.nii.gz');
outvol = fullfile(session_dir,d{runNum},[func '.anat.nii.gz']);
bbreg_out_file = fullfile(session_dir,d{runNum},[func '_bbreg.dat']); % registration file
[~,~] = system(['mri_vol2vol --mov ' invol ...
    ' --targ ' fullfile(SUBJECTS_DIR,subject_name,'mri/T1.mgz') ' --o ' ...
    outvol ' --reg ' bbreg_out_file ' --nearest']);
%% Remove temporary file
system(['rm ' fullfile(session_dir,d{runNum},'tmp.nii.gz')]);