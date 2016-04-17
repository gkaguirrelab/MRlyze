function create_nuisance_regressors(session_dir,runNum,r,physio,motion,parfile,contrasts)

% creates a nuisance regressor text file, based on physiological noise and
% motion
%
%   Usage:
%   create_nuisance_regressors(session_dir,physio,motion,feat_dir,parfile)
%
%   defaults:
%   physio = 1; % create physiological noise regressors
%   motion = 1; % create motion regressors
%   feat_dir = 'raw_f.feat';
%   parfile = 'prefiltered_func_data_mcf.par';
%
%   Written by Andrew S Bock Oct 2014

%% Setup initial variables
if ~exist('physio','var')
    physio = 1;
end
if ~exist('motion','var')
    motion = 1;
end
if ~exist('parfile','var')
    parfile = 'brf_motion_params.par';
end
%% Find bold run directories
d = listdir(fullfile(session_dir,'*BOLD_*'),'dirs');
if isempty(d)
    d = listdir(fullfile(session_dir,'*EPI_*'),'dirs');
end
if isempty(d)
    d = listdir(fullfile(session_dir,'RUN*'),'dirs');
end
disp(['Session_dir = ' session_dir]);
%disp(['Number of runs = ' num2str(nruns)]);
clear pulse newmotion
%% Create physiological noise regressors
if physio
    % Create physiological regressors
    outDir = fullfile(session_dir,d{r});
    dicomDir = fullfile(session_dir,'DICOMS',d{r});
    pulseDir = fullfile(session_dir,'PulseOx');
    pulsFile = fullfile(pulseDir,['RUN' num2str(runNum) '.puls']);
    pulse = PulseResp(dicomDir,pulsFile,outDir);
    % The above 'PulseResp' step created a file called 'puls.mat' in the
    % run directory
end
%% Create motion regressors, orthogonal to task contrasts
if motion
    % Create physiological regressors
    outDir = fullfile(session_dir,d{r});
    newmotion = regress_motion(fullfile(outDir,parfile),contrasts,outDir);
end
%% Write to file
outDir = fullfile(session_dir,d{r});
if physio && motion
    dlmwrite(fullfile(outDir,'nuisance_regressors.txt'),[pulse.all,newmotion],'delimiter','\t','precision','%10.5f')
elseif physio && ~motion
    dlmwrite(fullfile(outDir,'nuisance_regressors.txt'),pulse.all,'delimiter','\t','precision','%10.5f')
elseif ~physio && motion
    dlmwrite(fullfile(outDir,'nuisance_regressors.txt'),newmotion,'delimiter','\t','precision','%10.5f')
end
