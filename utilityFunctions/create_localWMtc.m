function [WMtc] = create_localWMtc(session_dir,runNum,func)
% Creates local white matter timecourses for each voxel, and saves these
%   timecourses as a 4D volume [func '.WMtc.nii.gz']. The default 'func'
%   input is 'brf', so the final 4D volume will be 'brf.WMtc.nii.gz'. Will
%   also return an output 4D matrix 'WMtc'.
%
% Usage:
%   [WMtc] = create_localWMtc(session_dir,runNum,func)
%
%   Written by Andrew S Bock Apr 2015

%% Set default parameters
if ~exist('session_dir','var')
    error('"session_dir" not defined')
end
if ~exist('func','var')
    func = 'brf'; % functional data file used for registration
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
SaveLogInfo(session_dir, mfilename,session_dir,runNum,func);

%% Create local WM timecourse
for rr = runNum;
    % Load in brain and white matter masks, as well as the functional timecourse file
    brain = load_nifti(fullfile(session_dir,d{rr},[func '.brainmask.nii.gz']));
    wmmask = load_nifti(fullfile(session_dir,d{rr},[func '.aseg.wm.nii.gz']));
    fmri = load_nifti(fullfile(session_dir,d{rr},[func '.tf.nii.gz']));
    voxsize = fmri.pixdim(2);
    dims=size(fmri.vol);
    tc = reshape(fmri.vol,dims(1)*dims(2)*dims(3),dims(4));
    tc = tc';
    % Compute the local white matter timecourses for each voxel
    disp('Finding localWM Voxels');
    localWMtc = local_WM(brain.vol,wmmask.vol,tc,voxsize);
    % Save the output volume
    disp('Saving local white matter volume');
    localWMtc = localWMtc';
    WMtc = reshape(localWMtc,size(fmri.vol));
    fmri.vol = WMtc;
    save_nifti(fmri,fullfile(session_dir,d{rr},[func '.WMtc.nii.gz']));
    disp('done.')
end