function register_functional_feat(check,adjust,session_dir,subject,func)
%   Registers functional runs in feat to freesurfer anatomical
%
%   Usage:
%   register_functional_feat(session_dir,subject,func,SUBJECTS_DIR)
%
%   e.g. register_functional_feat(1,0,'~/data/ASB'),'ASB')
%
%   Defaults:
%   check = 0; do not check registration
%   adjust = 0; do not manually adjust registration
%   session_dir = NO DEFAULT, must define
%   subject = NO DEFAULT, must define
%   func = 'filtered_func_data'
%   hemi = {'lh','rh'};
%
%   Outputs:
%   The function will display the minimum cost value of registraion, for
%   each feat direcotry
%
%   The function will also create the following files:
%   # = number, for each stat file (e.g. cope1, tstat1, zstat1)
%   m = hemisphere, either l - left, or r - right
%   1) mh_surf_cope#.nii.gz; cope, unsmoothed
%   2) mh_smooth_surf_cope#.nii.gz; cope, smoothed on surface, 5mm kernel
%   3) mh_surf_tstat#.nii.gz; tstat, unsmoothed
%   4) mh_smooth_surf_tstat#.nii.gz; tstat, smoothed on surface, 5mm kernel
%   5) mh_surf_zstat#.nii.gz; zstat, unsmoothed
%   6) mh_smooth_surf_zstat#.nii.gz; zstat, smoothed on surface, 5mm kernel
%
%   Written by Andrew S Bock Sept 2014
%% Set default variables
if ~exist('check','var')
    check = 0; % check registration
end
if ~exist('adjust','var')
    adjust = 0; % manually adjust registration
end
if ~exist('session_dir','var')
    error('"session_dir" not defined')
end
if ~exist('subject','var')
    error('"subject" not defined')
end
if ~exist('func','var')
    func = 'filtered_func_data'; % functional data file
end
if ~exist('hemi','var')
    hemi = {'lh','rh'};
end
if ~exist('stat_files','var')
    stat_files = {'cope','tstat','zstat'};
end
%% Set up FSL variables
fsl_path = '/usr/local/fsl/';
setenv('FSLDIR',fsl_path)
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
curpath = getenv('PATH');
setenv('PATH',sprintf('%s:%s',fullfile(fsl_path,'bin'),curpath));
%% Find feat directories in the session directory
d = listdir(fullfile(session_dir,'*.feat'),'dirs');
nruns = length(d);
disp(['found ' num2str(nruns) ' feat directories']);
%% Register functional to freesurfer anatomical
progBar = ProgressBar(nruns, 'Registering functional runs to freesurfer anatomical...');
mincost = nan(nruns,1);
for r = 1:nruns
    filefor_reg = fullfile(session_dir,d{r},[func '.nii.gz']); % Functional file for bbregister
    bbreg_out_file = fullfile(session_dir,d{r},'bbreg.dat'); % name registration file
    % Register functional volume to anatomical volume using bbregister
    [~,~] = system(['bbregister --s ' subject ' --mov ' filefor_reg ' --reg ' bbreg_out_file ' --init-fsl --t2']);
    if check % Check the registration
        [~,~] = system(['tkregister2 --mov ' filefor_reg ' --reg ' bbreg_out_file ' --surf']);
    end
    if adjust % reiterate following manual adjustment
        [~,~] = system(['bbregister --s ' subject ' --mov ' filefor_reg ' --reg ./' ...
            bbreg_out_file ' --init-reg ' bbreg_out_file ' --t2']);
        [~,~] = system(['tkregister2 --mov ' filefor_reg ' --reg ./' bbreg_out_file ' --surf']);
    end
    progBar(r);
    load(fullfile(session_dir,d{r},'bbreg.dat.mincost'));
    mincost(r) = bbreg_dat(1);
    clear bbreg_dat
end
%% Display the results of registration
disp('The min cost values of registration for each feat directory are:')
disp(num2str(mincost));
%% Project cope files to surface
% use the bbreg_out_file registration file created above
progBar = ProgressBar(nruns, 'Projecting stat files to surface...');
for r = 1:nruns
    % For each stat file (e.g. cope, tstat, zstat)
    for s = 1:length(stat_files)
        % Find the number of stat files in each feat directory
        f = listdir(fullfile(session_dir,d{r},'stats',[stat_files{s} '*']),'files');
        nfiles = length(f);
        for n = 1:nfiles
            cd(fullfile(session_dir,d{r},'stats'));
            % Define registration file
            bbreg_out_file = fullfile(session_dir,d{r},'bbreg.dat');
            % Project cope files to surface, for each hemisphere
            for h = 1:length(hemi)
                % Unsmoothed - native
                [~,~] = system(['mri_vol2surf --src ./' f{n} ' --reg ' bbreg_out_file ' --hemi ' hemi{h} ...
                    ' --out ./' hemi{h} '_surf_' f{n} ' --projfrac 0.5']);
                % Unsmooth - fsaverage_sym
                [~,~] = system(['mri_surf2surf --hemi ' hemi{h} ' --srcsubject ' subject ' --srcsurfval ' ...
                    './' hemi{h} '_surf_' f{n} ' --trgsubject fsaverage_sym --trgsurfval ' ...
                    './' hemi{h} '_fsavg_symsurf_' f{n}]);
                % Smoothed - native
                [~,~] = system(['mri_vol2surf --src ./' f{n} ' --reg ' bbreg_out_file ' --hemi ' hemi{h} ...
                    ' --surf-fwhm 5 --out ./' hemi{h} '_smooth_surf_' f{n} ' --projfrac 0.5']);
                % Smoothed - fsaverage_sym
                [~,~] = system(['mri_surf2surf --hemi ' hemi{h} ' --srcsubject ' subject ' --srcsurfval ' ...
                    './' hemi{h} '_smooth_surf_' f{n} ' --trgsubject fsaverage_sym --trgsurfval ' ...
                    './' hemi{h} '_smooth_fsavg_symsurf_' f{n}]);
            end
        end
    end
    progBar(r);
end
