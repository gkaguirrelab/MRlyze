function projectCopes2anat(session_dir,subject_name,runNums,func)

% Projects 'cope*.nii.gz' files in a FEAT directory to anatomical space
%
%   Usage:
%   projectCopes2anat(session_dir,subject_name,runNums,func);
%
%   Written by Andrew S Bock May 2016

%% set defaults
if ~exist('func','var')
    func = 'wdrf.tf';
end
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
%% Find bold directories
b = find_bold(session_dir);

%% Project cope files to anatomical space
for i = runNums
    % Get directory and cope files
    inDir = fullfile(session_dir,b{i},[func '.feat'],'stats');
    copeFiles = listdir(fullfile(inDir,'cope*.nii.gz'),'files');
    targ_vol = fullfile(SUBJECTS_DIR,subject_name,'mri','T1.mgz');
    % Project mean functional volume to anatomical space
    meanfile = fullfile(session_dir,b{i},[func '.feat'],'mean_func.nii.gz');
    meanout = fullfile(session_dir,b{i},[func '.feat'],'mean_func.anat.nii.gz');
    tmpreg = listdir(fullfile(session_dir,b{i},'*bbreg.dat'),'files');
    bbreg_out_file = fullfile(session_dir,b{i},tmpreg{1}); % name registration file
    system(['mri_vol2vol --mov ' meanfile ...
        ' --targ ' targ_vol ' --reg ' bbreg_out_file ...
        ' --o ' meanout ' --trilin']);
    % Project cope files to anatomical space
    for j = 1:length(copeFiles)
        if isempty(strfind(copeFiles{j},'anat'))
            % Project cope to anatomical space
            copefile = fullfile(inDir,copeFiles{j});
            extInd = strfind(copeFiles{j},'.nii.gz');
            copeout = fullfile(inDir,[copeFiles{j}(1:extInd) 'anat.nii.gz']);
            tmpreg = listdir(fullfile(session_dir,b{i},'*bbreg.dat'),'files');
            bbreg_out_file = fullfile(session_dir,b{i},tmpreg{1}); % name registration file
            system(['mri_vol2vol --mov ' copefile ...
                ' --targ ' targ_vol ' --reg ' bbreg_out_file ...
                ' --o ' copeout ' --trilin']);
        end
    end
end