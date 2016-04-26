function run_pRF(session_dir,subject_name,runNum,hemi,srcROI,func,fieldSize,imFileName,paramsFileName)

% Run the population receptive field (pRF) pipeline
%
%   Usage:
%   run_pRF(session_dir,subject_name,runNum,hemi,srcROI)
%
%   Written by Andrew S Bock Jun 2015

%% Find bold run directories
d = find_bold(session_dir);
%% Set up defaults
if ~exist('runNum','var')
    runNum = 1:length(d);
end
if ~exist('hemi','var')
    hemi = 'lh';
end
if ~exist('srcROI','var')
    srcROI = 'volume';
end
if ~exist('fieldSize','var')
    fieldSize = 6.2116; % UPenn SC7T; 10.4721 - UPenn SC3T
end
if ~exist('imFileName','var')
    imFileName = 'bars_images.mat';
end
if ~exist('paramsFileName','var');
    paramsFileName = 'bars_params.mat';
end
%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,subject_name,runNum,hemi,srcROI,func,fieldSize,imFileName,paramsFileName)

%% Run pRF fitting on surface
for rr = runNum;
    switch srcROI
        case 'cortex'
            srcfile = fullfile(session_dir,d{rr},[func '.' hemi '.nii.gz']);
            areas = load_nifti(fullfile(session_dir,'pRFs',...
                'anat_templates',[hemi '.areas.nii.gz']));
            srcind = 1:length(areas.vol); % entire cortex
        case 'volume'
            srcfile = fullfile(session_dir,d{rr},[func '.nii.gz']);
            %binfile = fullfile(session_dir,d{rr},'brf.aseg.gm.nii.gz');
            binfile = fullfile(session_dir,d{rr},'single_TR.nii.gz');
            src = load_nifti(binfile);
            srcind = find(src.vol > 0 & ~isnan(src.vol));
        case 'LGN'
            srcfile = fullfile(session_dir,d{rr},[func '.nii.gz']);
            binfile = fullfile(session_dir,d{rr},[hemi '.LGN.nii.gz']);
            src = load_nifti(binfile);
            srcind = find(src.vol > 0 & ~isnan(src.vol));
    end
    imFile = fullfile(session_dir,'Stimuli',['run' num2str(rr)],imFileName);
    paramsFile = fullfile(session_dir,'Stimuli',['run' num2str(rr)],paramsFileName);
    % Run the CF pipeline
    disp(['Generating ' srcROI ' pRFs for ' d{rr}]);
    do_pRF(session_dir,subject_name,rr,hemi,srcROI,func,srcind,srcfile,fieldSize,imFile,paramsFile)
end
