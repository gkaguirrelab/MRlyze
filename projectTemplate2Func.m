function [roiInd] = projectTemplate2Func(params)

% Projects an anatomical template to functional space, and output the 
%   indices of the specified ROI
%
%   Usage:
%       [roiInd] = projectTemplate2Func(params)
%
%   Required:
%       params.sessionDir   = '/path/to/session/directory'
%
%   Defaults:
%       params.runNum       = 1; % first bold directory
%       params.roiType      = 'V1'; % could also be 'V2V3' and 'LGN'
%       params.func         = 'wdrf.tf';
%       params.eccRange     = [2.5 32]; % based on MaxMel data
%
%   Written by Andrew S Bock Oct 2016

%% Set defaults
if ~isfield(params,'runNum');
    params.runNum           = 1;
end
if ~isfield(params,'roiType');
    params.roiType          = 'V1';
end
if ~isfield(params,'func');
    params.func             = 'wdrf.tf';
end
if ~isfield(params,'eccRange');
    params.eccRange         = [2.5 32]; % based on MaxMel data
end
bbregName                   = 'func_bbreg.dat';
boldDirs                    = find_bold(params.sessionDir);
%% Set the file names
switch params.roiType
    case {'V1' 'V2V3'}
        areaInName          = 'mh.areas.anat.vol.nii.gz';
        areaOutName         = 'mh.areas.func.vol.nii.gz';
        eccInName           = 'mh.ecc.anat.vol.nii.gz';
        eccOutName          = 'mh.ecc.func.vol.nii.gz';
    case 'LGN'
        areaInName          = 'mh.LGN.nii.gz';
        areaOutName         = 'mh.LGN.func.vol.nii.gz';
end
%% Project anatomical template to functional space
% Project area anat to func
areaInFile                  = fullfile(params.sessionDir,'anat_templates',areaInName);
bbregFile                   = fullfile(params.sessionDir,boldDirs{params.runNum},bbregName); % registration file
areaOutFile                 = fullfile(params.sessionDir,boldDirs{params.runNum},areaOutName);
system(['mri_vol2vol --mov ' fullfile(params.sessionDir,boldDirs{params.runNum},[params.func '.nii.gz']) ...
    ' --targ ' areaInFile ' --o ' areaOutFile ...
    ' --reg ' bbregFile ' --inv --nearest']);
areaData                    = load_nifti(areaOutFile);
% If 'V1stim' or 'V2V3stim', also project eccentricity file
switch params.roiType
    case {'V1' 'V2V3'}
        eccInFile           = fullfile(params.sessionDir,'anat_templates',eccInName);
        bbregFile           = fullfile(params.sessionDir,boldDirs{params.runNum},bbregName); % registration file
        eccOutFile          = fullfile(params.sessionDir,boldDirs{params.runNum},eccOutName);
        system(['mri_vol2vol --mov ' fullfile(params.sessionDir,boldDirs{params.runNum},[params.func '.nii.gz']) ...
            ' --targ ' eccInFile ' --o ' eccOutFile ...
            ' --reg ' bbregFile ' --inv --nearest']);
        eccData             = load_nifti(eccOutFile);
end
switch params.roiType
    case 'LGN'
        roiInd              = find(abs(areaData.vol)==1);
    case 'V1'
        roiInd              = find(abs(areaData.vol)==1 & ...
            eccData.vol>params.eccRange(1) & eccData.vol<params.eccRange(2));
    case 'V2V3'
        roiInd              = find((abs(areaData.vol)==2 | abs(areaData.vol)==3) & ...
            eccData.vol>params.eccRange(1) & eccData.vol<params.eccRange(2));
end