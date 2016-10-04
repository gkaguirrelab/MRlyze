function [roiInd] = getROI(params)

% Outputs the ROI indices specifed by `params` inputs
%
%   Usage:
%       [roiInd] = getROI(params)
%
%   Required:
%       params.sessionDir   = '/path/to/session/directory'
%
%   Defaults:
%       params.runNum       = 1; % first bold directory
%       params.roiType      = 'V1'; % could also be 'V2V3' and 'LGN'
%       params.eccRange     = [0 30]; % based on MaxMel data
%
%   Written by Andrew S Bock Oct 2016

%% Set defaults
if ~isfield(params,'runNum');
    params.runNum           = 1;
end
if ~isfield(params,'roiType');
    params.roiType          = 'V1';
end
if ~isfield(params,'eccRange');
    params.eccRange         = [0 30]; % based on MaxMel data
end
boldDirs                    = find_bold(params.sessionDir);
%% Set the file names
switch params.roiType
    case {'V1' 'V2V3'}
        areaName            = 'mh.areas.func.vol.nii.gz';
        eccName             = 'mh.ecc.func.vol.nii.gz';
    case 'LGN'
        areaName            = 'mh.LGN.func.vol.nii.gz';
end
%% Load files
areaFile                    = fullfile(params.sessionDir,boldDirs{params.runNum},areaName);
eccFile                     = fullfile(params.sessionDir,boldDirs{params.runNum},eccName);
areaData                    = load_nifti(areaFile);
eccData                     = load_nifti(eccFile);
%% Get ROI indices
roiInd                      = find(abs(areaData.vol)==1 & ...
                            eccData.vol>params.eccRange(1) & ...
                            eccData.vol<params.eccRange(2));