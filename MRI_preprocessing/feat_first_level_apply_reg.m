function outdir = feat_first_level_apply_reg(session_dirs,subject_name,dirs,func,BACKUP,OVERWRITE_REG_STANDARD)

%   Following feat_stats, this function will apply registrations so that
%   higher level feat will run.
%
%   inputs:
%   session_dirs - cell containing strings of the paths to the session
%   directories of interest
%   subject_name - freesurfer subject name
%   dirs - vector of feat directories (e.g. [1 4 7 10 13 16])
%   directories of interest
%   func - functional data file (default - 'dbrf.tf')
%
%   Written by Andrew S Bock Dec 2014
%
%   3/4/15  ms      Pulled out of feat_higher_level

%% Set default parameters
if ~exist('session_dirs','var')
    error('"session_dirs" not defined')
end
if ~exist('subject_name','var')
    error('"subject_name" not defined')
end
if ~exist('dirs','var')
    error('"dirs" not defined') % feat directories
end
if ~exist('func','var')
    func = 'dbrf.tf'; % functional data file
end
if ~exist('BACKUP', 'var')
    BACKUP = false;
end
if ~exist('OVERWRITE_REG_STANDARD', 'var')
    OVERWRITE_REG_STANDARD = false;
end

%% Loop through session directories
for s = 1:length(session_dirs)
    session_dir = session_dirs{s};
    % Find bold run directories
    d = listdir(fullfile(session_dir,'*BOLD_*'),'dirs');
    if isempty(d)
        d = listdir(fullfile(session_dir,'*EPI_*'),'dirs');
    end
    if isempty(d)
        d = listdir(fullfile(session_dir,'RUN*'),'dirs');
    end
    nruns = length(d);
    % Copy over bbregister registration file
    % overwrite the example_func2standard.mat in each feat dir
    for r = 1:nruns
        ct = ct+1;
        feat_dirs{ct} = fullfile(session_dir,d{r},feat_dir);
        % Backup the FSL registration (which is worse than bbregister)
        if BACKUP
                copyfile(fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.mat'),...
                    fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.mat_BAK'));
                copyfile(fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.nii.gz'),...
                    fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.nii.gz_BAK'));
                copyfile(fullfile(session_dir,d{r},feat_dir,'reg','standard.nii.gz'),...
                    fullfile(session_dir,d{r},feat_dir,'reg','standard.nii.gz_BAK'));
                copyfile(fullfile(session_dir,d{r},feat_dir,'reg','standard2example_func.mat'),...
                    fullfile(session_dir,d{r},feat_dir,'reg','standard2example_func.mat_BAK'));
                copyfile(fullfile(session_dir,d{r},feat_dir,'reg_standard','example_func.nii.gz'),...
                    fullfile(session_dir,d{r},feat_dir,'reg_standard','example_func.nii.gz_BAK'));
                copyfile(fullfile(session_dir,d{r},feat_dir,'reg_standard','mask.nii.gz'),...
                    fullfile(session_dir,d{r},feat_dir,'reg_standard','mask.nii.gz_BAK'));
                copyfile(fullfile(session_dir,d{r},feat_dir,'reg_standard','mean_func.nii.gz'),...
                    fullfile(session_dir,d{r},feat_dir,'reg_standard','mean_func.nii.gz_BAK'));
        end
        % Overwrite FSL registration with Freesurfer bbregister registration
        system(['tkregister2 --mov ' fullfile(session_dir,d{r},[func '.nii.gz']) ...
            ' --reg ' fullfile(session_dir,d{r},'brf_bbreg.dat') ' --fslregout ' ...
            fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.mat') ...
            ' --surfs --noedit']);
        % Overwrite inverse registration matrix
        system(['convert_xfm -omat ' fullfile(session_dir,d{r},feat_dir,'reg','standard2example_func.mat') ...
            ' -inverse ' fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.mat')]);
        % Overwrite standard
        system(['mri_convert ' fullfile(SUBJECTS_DIR,subject_name,'mri','brain.mgz') ' ' ...
            fullfile(session_dir,d{r},feat_dir,'reg','standard.nii.gz')]);
        % Overwrite example_func2standard.nii.gz
        system(['flirt -in ' fullfile(session_dir,d{r},feat_dir,'reg','example_func.nii.gz') ...
            ' -ref ' fullfile(session_dir,d{r},feat_dir,'reg','standard.nii.gz') ' -out ' ...
            fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.nii.gz') ...
            ' -init ' fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.mat') ...
            ' -applyxfm']);
        % Overwrite reg_standard files
        if OVERWRITE_REG_STANDARD
                copyfile(fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.nii.gz'),...
                    fullfile(session_dir,d{r},feat_dir,'reg_standard','example_func.nii.gz'));
                copyfile(fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.nii.gz'),...
                    fullfile(session_dir,d{r},feat_dir,'reg_standard','mean_func.nii.gz'));
                system(['fslmaths ' fullfile(session_dir,d{r},feat_dir,'reg_standard','example_func.nii.gz') ...
                    ' -bin ' fullfile(session_dir,d{r},feat_dir,'reg_standard','mask.nii.gz')]);
        end
    end
    % Save 1st session, first 'dirs' directory for later
    if s == 1
        outdir = fullfile(session_dir,d{dirs(1)});
    end
end