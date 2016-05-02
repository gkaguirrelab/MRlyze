function overwrite_feat_reg(session_dir,subject_name,runNums,featName,SUBJECTS_DIR)

% Overwrites the feat registration matrices, as well as the resulting 
%   registration volumes, using the bbregister registration created using
%   'register_feat'
%
%   Usage:
%   overwrite_feat_reg(session_dir,subject_name,runNums,featName,SUBJECTS_DIR)
%
%   Written by Andrew S Bock Mar 2015

%% Set default parameters
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
feat_dir = [featName '.feat'];
%% Loop through session directories
% Find bold run directories
d = find_bold(session_dir);
% Copy over bbregister registration file
% overwrite the example_func2standard.mat in each feat dir
for rr = runNums
    % Overwrite FSL registration with Freesurfer bbregister registration
    tmpreg = listdir(fullfile(session_dir,d{rr},'*bbreg.dat'),'files');
    fileforreg = fullfile(session_dir,d{rr},tmpreg{1});
    system(['tkregister2 --mov ' fullfile(session_dir,d{rr},'rf.nii.gz') ...
        ' --reg ' fileforreg ' --fslregout ' ...
        fullfile(session_dir,d{rr},feat_dir,'reg','example_func2standard.mat') ...
        ' --surfs --noedit']);
    % Overwrite inverse registration matrix
    system(['convert_xfm -omat ' fullfile(session_dir,d{rr},feat_dir,'reg','standard2example_func.mat') ...
        ' -inverse ' fullfile(session_dir,d{rr},feat_dir,'reg','example_func2standard.mat')]);
    % Overwrite standard
    system(['mri_convert ' fullfile(SUBJECTS_DIR,subject_name,'mri','brain.mgz') ' ' ...
        fullfile(session_dir,d{rr},feat_dir,'reg','standard.nii.gz')]);
    % Overwrite example_func2standard.nii.gz
    system(['flirt -in ' fullfile(session_dir,d{rr},feat_dir,'reg','example_func.nii.gz') ...
        ' -ref ' fullfile(session_dir,d{rr},feat_dir,'reg','standard.nii.gz') ' -out ' ...
        fullfile(session_dir,d{rr},feat_dir,'reg','example_func2standard.nii.gz') ...
        ' -init ' fullfile(session_dir,d{rr},feat_dir,'reg','example_func2standard.mat') ...
        ' -applyxfm']);
end