function remove_noise_fast(session_dir,subject_name,runNum,func,remove_task,anat,motion,physio,SUBJECTS_DIR)
% Removes physiological and other non-neuronal noise.
%
%   Usage:
%   remove_noise_fast(session_dir,subject_name,runNum,func,remove_task,anat,motion,physio,SUBJECTS_DIR)
%
%   Written by Andrew S Bock Apr 2015

%% Set default parameters
if ~exist('session_dir','var')
    error('"session_dir" not defined')
end
if ~exist('subject_name','var')
    error('"subject_name" not defined')
end
if ~exist('func','var')
    func = 'brf'; % functional data file used for registration
end
if ~exist('remove_task','var')
    remove_task = 0; % if '1', regresses out task conditions from motion
end
if ~exist('anat','var')
    anat = 1;  % if '0', won't remove anatomical ROI signals
end
if ~exist('motion','var')
    motion = 1;  % if '0', won't remove motion (and pulse) signals
end
if ~exist('physio','var')
    physio = 1;  % if '0', won't remove motion (and pulse) signals
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
%% In case no noise removal to be done, just return
if ~remove_task && ~anat && ~motion && ~physio
    return
end
%% Find bold run directories
d = find_bold(session_dir);
nruns = length(d);
disp(['Session_dir = ' session_dir]);
disp(['Number of runs = ' num2str(nruns)]);
%% Set runs
if ~exist('runNum','var')
    runNum = 1:length(d);
end
%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,subject_name,runNum,func,remove_task,anat,motion,SUBJECTS_DIR);

%% Load in timecourses
for rr = runNum
    disp(['Load functional timecoures and ROIs from ' fullfile(session_dir,d{rr}) '...']);
    % navigate to the run directory
    cd(fullfile(session_dir,d{rr}));
    % Load the timecourse for this run, already temporally filtered
    fmri = load_nifti(fullfile(session_dir,d{rr},[func '.tf.nii.gz']));
    dims=size(fmri.vol);
    tc = reshape(fmri.vol,dims(1)*dims(2)*dims(3),dims(4));
    tc = tc';
    % Load Brain mask
    brain = load_nifti(fullfile(session_dir,d{rr},[func '.brainmask.nii.gz']));
    %% Load anatomical regressors
    if ~anat
        anat_noise = [];
    else
        anat_noise = load(fullfile(session_dir,d{rr},'anat_params.txt'));
    end
    %% Load up other nuisance regressors
    if ~motion
        motion_noise = [];
    else
        motion_noise = load(fullfile(session_dir,d{rr},'motion_params.txt'));
        if remove_task
            % Remove task
            bold_dir = fullfile(session_dir,d{rr});
            TR = fmri.pixdim(5)/1000;
            if TR < 0.1
                error('TR is less than 0.1, most likely input nifti TR not in msec')
            end
            lengthTC = size(tc,1);
            % Convert task conditions in to timecoures (convolve with HRF) that
            % are at the resolution of the TR.
            [outTC] = convert_task2tc(bold_dir,TR,lengthTC);
            % regress out task from motion
            motion_noise = regress_task(motion_noise,outTC);
        end
    end
    if ~physio
        physio_noise = [];
    else
        physio_noise = load(fullfile(session_dir,d{rr},'pulse_params.txt'));
    end
    noise = [anat_noise,motion_noise,physio_noise];
    % remove any means and linear trends
    noise = detrend(noise);
    %% Calculate global signal
    GB.ind = find(brain.vol);
    GB.tc = mean(tc(:,GB.ind),2);
    GB.tc = detrend(GB.tc); % remove mean and linear trend
    %% Regress out noise
    newtc = tc;
    Gnewtc = tc;
    varexp.all = zeros(length(tc),1);
    varexp.Gall = zeros(length(tc),1);
    regressMat = noise;
    [orthmat] = orth(regressMat);
    [Gorthmat] = orth([regressMat,GB.tc]);
    % Remove noise, and noise+global signal
    [oB] = [ones(size(tc,1),1),orthmat]\tc;
    [GoB] = [ones(size(tc,1),1),Gorthmat]\tc;
    orthbeta = oB(2:end,:);
    Gorthbeta = GoB(2:end,:);
    newtc = tc-orthmat*(orthbeta);
    Gnewtc = tc-Gorthmat*(Gorthbeta);
    disp('done.');
    % Save timecourses
    disp('Saving filtered timecourses...');
    fmri = load_nifti(fullfile(session_dir,d{rr},[func '.tf.nii.gz']));
    newtc = newtc';
    Gnewtc = Gnewtc';
    newtc = reshape(newtc,size(fmri.vol));
    Gnewtc = reshape(Gnewtc,size(fmri.vol));
    % load template file, for use in saving outputs
    fmri.vol = newtc;
    dsavename = ['d' func '.tf.nii.gz'];
    save_nifti(fmri,fullfile(session_dir,d{rr},dsavename));
    gdsavename = ['gd' func '.tf.nii.gz'];
    fmri.vol = Gnewtc;
    save_nifti(fmri,fullfile(session_dir,d{rr},gdsavename));
    disp('done.');
end