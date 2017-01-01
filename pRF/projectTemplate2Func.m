function projectTemplate2Func(params)

% Projects anatomical template(s) to functional space
%
%   Usage:
%       projectTemplate2Func(params)
%
%   Required:
%       params.sessionDir   = '/path/to/session/directory'
%
%   Defaults:
%       params.runNum       = 1; % first bold directory
%       params.roiType      = 'V1'; % could also be 'V2V3' and 'LGN'
%       params.func         = 'wdrf.tf';
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
bbregName                   = 'func_bbreg.dat';
boldDirs                    = find_bold(params.sessionDir);
%% Set the file names
switch params.roiType
    case {'V1' 'V2V3'}
        areaInName          = 'mh.areas.anat.vol.nii.gz';
        areaOutName         = 'mh.areas.func.vol.nii.gz';
        eccInName           = 'mh.ecc.anat.vol.nii.gz';
        eccOutName          = 'mh.ecc.func.vol.nii.gz';
        polInName           = 'mh.pol.anat.vol.nii.gz';
        polOutName          = 'mh.pol.func.vol.nii.gz';
    case 'LGN'
        areaInName          = 'mh.LGN.nii.gz';
        areaOutName         = 'mh.LGN.func.vol.nii.gz';
end
%% Project area template to functional space
areaInFile                  = fullfile(params.sessionDir,'anat_templates',areaInName);
bbregFile                   = fullfile(params.sessionDir,boldDirs{params.runNum},bbregName); % registration file
areaOutFile                 = fullfile(params.sessionDir,boldDirs{params.runNum},areaOutName);
system(['mri_vol2vol --mov ' fullfile(params.sessionDir,boldDirs{params.runNum},[params.func '.nii.gz']) ...
    ' --targ ' areaInFile ' --o ' areaOutFile ...
    ' --reg ' bbregFile ' --inv --nearest']);
%% If 'V1' or 'V2V3', also project eccentricity and polar angle files
switch params.roiType
    case {'V1' 'V2V3'}
        % eccentricity
        eccInFile           = fullfile(params.sessionDir,'anat_templates',eccInName);
        bbregFile           = fullfile(params.sessionDir,boldDirs{params.runNum},bbregName); % registration file
        eccOutFile          = fullfile(params.sessionDir,boldDirs{params.runNum},eccOutName);
        system(['mri_vol2vol --mov ' fullfile(params.sessionDir,boldDirs{params.runNum},[params.func '.nii.gz']) ...
            ' --targ ' eccInFile ' --o ' eccOutFile ...
            ' --reg ' bbregFile ' --inv --nearest']);
        % polar angle
        polInFile           = fullfile(params.sessionDir,'anat_templates',polInName);
        bbregFile           = fullfile(params.sessionDir,boldDirs{params.runNum},bbregName); % registration file
        polOutFile          = fullfile(params.sessionDir,boldDirs{params.runNum},polOutName);
        system(['mri_vol2vol --mov ' fullfile(params.sessionDir,boldDirs{params.runNum},[params.func '.nii.gz']) ...
            ' --targ ' polInFile ' --o ' polOutFile ...
            ' --reg ' bbregFile ' --inv --nearest']);
end