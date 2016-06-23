%% Template script for analyzing MRI data using MRlyze
%
%   Required software:
%       Freesurfer, FSL
%       IMPORTANT:
%           start matlab from terminal, so environmental variables,
%               libraries are set correctly
%           add $FREESURFER_HOME/matlab to your matlab path
%
%   Directory structure:
%       data_directory -> project directory -> subject directory ->
%           session directory -> dicom directory
%
%           e.g. ~/data/Retinotopy/ASB/10012014/DICOMS
%
%       If physiological measures are collected, using the following
%           directory structure:
%
%       data_directory -> project directory -> subject directory ->
%           session directory -> physio directory
%
%           e.g. ~/data/Retinotopy/ASB/10012014/PulseOx
%
%   I recommend creating a project specific master file in a separate
%   directory. For example, copy this MASTER.m file to a project
%   specific folder:
%
%   /User/Shared/Matlab/<project_name>/<project_name>_MASTER.m
%
%   Written by Andrew S Bock Jun 2016

%% Define global paths
% Change these paths according to the environment currently in use.
% output_dir is usually the link to the lab's dropbox. data_dir must have
% the standard initial data organization. 'matDir' is the parent
% directory where stimulus .mat files are saved during the experiment.
dataDir                 = '/data/jag/MELA/MOUNT_SINAI/';
projectGitDir           = '/data/jag/MELA/Matlab/gkaguirrelab/SC-7T/';
dropboxDir              = getenv('DROPBOX_DIR');
outputDropboxDir        = fullfile(dropboxDir,'MELA_analysis','HCLV_Photo_7T');
projectName             = 'MELA_data'; % Project data name
olprotocolName          = 'HCLV_Photo'; % Protocol Name
scannerProtocolName     = 'HCLV_Photo_7T'; % Protocol Name
stimDir                 = fullfile(dropboxDir,projectName,scannerProtocolName);
logDir                  = getenv('LOG_DIR'); % log files on the UPenn cluster
SUBJECTS_DIR            = getenv('SUBJECTS_DIR'); % freesurfer subject directory
saveFigDir              = fullfile(outputDropboxDir,'MOUNT_SINAI_TTFs');
%% Define subject variables
% Subject names
subjects = {...
    'subject1' ...
    'subject2' ...
    };
% Session dates
session_dates = {...
    '041416' ...
    '041516' ...
    };
% Freesurfer subject names
fsSubNames = {...
    'HERO_subject1_7T' ...
    'HERO_subject2_7T' ...
    };
% number of bold runs in each session directory
allRuns = [12 26];
% One Light Protocol names
olprotocolNames = {...
    'MelanopsinMRMaxMel' ...
    'MelanopsinMRMaxMel' ...
    };
% define 'wrap around' runs for each session
wrapAround = cell(1,length(allRuns));
for i = 1:length(allRuns)
    wrapAround{i} = zeros(1,allRuns(i));
end
wrapAround{1}([1, 2, 3, 7, 8, 9]) = 1;
wrapAround{2}([1, 2, 3, 7, 8, 9, 13, 14, 15, 19, 20, 21, 26]) = 1;
% define any 'excluded' runs (e.g. poor task performance), if any
excludeRuns = {...
    [0 0 0 0 0 1 0 0 1 0 0 0] ... % two runs are exluded due to task performance
    zeros(1,allRuns(2)) ...
    };
% Functional volumes
funcs = {...
    'wdrf.tf' ...
    's5.wdrf.tf' ...
    };
% Stimulus Conditions
conditionNames = {...
    '_LightFlux_' ...
    '_L_minus_M_' ...
    '_S_'...
    };
plotNames = {'LightFlux' 'LminusM' 'S'};
% Stimulus Frequencies
freqs = {'2Hz' '4Hz' '8Hz' '16Hz' '32Hz' '64Hz'};
% ROIs
ROIs = {'V1all' 'V1low' 'V1mid' 'V1high' 'V2andV3' 'LGN' 'SC'};
% Hemispheres
hemis = {'lh' 'rh' 'mh'};
% Cope names
copeNames = {...
    'Hz2' ...
    'Hz4' ...
    'Hz8' ...
    'Hz16' ...
    'Hz32' ...
    'Hz64' ...
    };
copeNums = 1:length(copeNames);
% Cope threshold
pthreshVal = 0.05;
%% setup directories
sessions = cell(1,length(subjects));
stimFileDirs = cell(1,length(subjects));
protNames = cell(1,length(subjects));
stimFileNames = cell(1,length(subjects));
outStimuli = cell(1,length(subjects));
outPerformance = cell(1,length(subjects));



outFirstLevelStats = cell(1,length(subjects));
outF2p = cell(1,length(subjects));
outFishersChiSquared = cell(1,length(subjects));
outpsc_cope = cell(1,length(subjects));



for ss = 1:length(subjects)
    % session directory
    sessions{ss} = fullfile(dataDir,subjects{ss},session_dates{ss});
    % Stimulus file directory
    stimFileDirs{ss} = fullfile(stimDir,subjects{ss},session_dates{ss},'MatFiles');
    protNames{ss} = [subjects{ss} '-' scannerProtocolName];
    for j = 1:allRuns(ss)
        stimFileNames{ss}{j} = sprintf([protNames{ss} '-%02d.mat'],j);
    end
    % Stimulus files
    outStimuli{ss} = fullfile(sessions{ss},'Stimuli');
    if ~isdir(outStimuli{ss})
        mkdir(outStimuli{ss});
    end
    % Subject task performance
    outPerformance{ss} = fullfile(sessions{ss},'Performance');
    if ~isdir(outPerformance{ss})
        mkdir(outPerformance{ss});
    end
end
%% Generate regressors for statistical analysis
%   Also check subject task performance
for ss = 1:length(sessions)
    clear hits hTotal falseAlarms fTotal
    for mm = 1:allRuns(ss);
        % Get the stimulus files for each run
        matFile = fullfile(stimFileDirs{ss},stimFileNames{ss}{mm});
        % Create the text file regressors
        ol_regressors(matFile,outStimuli{ss},olprotocolName,wrapAround{ss}(mm));
        % Check performance
        [hits(mm),hTotal(mm),falseAlarms(mm),fTotal(mm)] = ...
            check_performance(matFile,olprotocolName);
    end
    % Save performance values
    fid = fopen(fullfile(outPerformance{ss}, 'performance.csv'), 'w');
    fprintf(fid, 'Hits,N,False alarms,N\n');
    fclose(fid);
    dlmwrite(fullfile(outPerformance{ss}, 'performance.csv'), ...
        [hits' hTotal' falseAlarms' fTotal'], '-append');
    fid = fopen(fullfile(outPerformance{ss}, 'performance.csv'), 'a');
    fprintf(fid, 'Total\n');
    fclose(fid);
    dlmwrite(fullfile(outPerformance{ss}, 'performance.csv'), ...
        [sum(hits) sum(hTotal) sum(falseAlarms) sum(fTotal)], '-append');
    fid = fopen(fullfile(outPerformance{ss}, 'performance.csv'), 'a');
    fprintf(fid, 'Percentages\n');
    fprintf(fid, '%.3f,,%.3f,', 100*(sum(hits)/sum(hTotal)), ...
        100*(sum(falseAlarms)/sum(fTotal)));
    fclose(fid);
end
%% Run stats
for ss = 1:length(sessions)
    d = find_bold(sessions{ss});
    for j = 1:length(d)
        for ff = 1:length(funcs)
            % load functional volume
            funcVol = fullfile(sessions{ss},d{j},[funcs{ff} '.nii.gz']);
            tmpV = load_nifti(funcVol);
            allTCs = reshape(tmpV.vol,size(tmpV.vol,1)*size(tmpV.vol,2)*size(tmpV.vol,3),size(tmpV.vol,4));
            TR = tmpV.pixdim(5); % TR in msec
            % Get the stimulus timing files to be used as regressors
            stimuli_dirs = listdir(outStimuli{ss},'dirs');
            tmpFiles = listdir(fullfile(outStimuli{ss},stimuli_dirs{j},'*_valid.txt'),'files');
            % Pull out the regressor files for each frequency
            for fq = 1:length(freqs)
                freqName = freqs{fq};
                for i = 1:length(tmpFiles)
                    tmp = strfind(tmpFiles{i},['_' freqName '_']);
                    if ~isempty(tmp)
                        EVs{fq} = fullfile(outStimuli{ss},stimuli_dirs{j},tmpFiles{i});
                    end
                end
            end
            % Add the attention task regressor
            taskFile = listdir(fullfile(outStimuli{ss},stimuli_dirs{j},'*attentionTask.txt'),'files');
            EVs{length(EVs)+1} = fullfile(outStimuli{ss},stimuli_dirs{j},taskFile{1});
            % Determine if run started with a 'wraparound' block
            WA = wrapAround{ss}(j);
            if WA
                WAFile = listdir(fullfile(outStimuli{ss},stimuli_dirs{j},'*wrapAround.txt'),'files');
                EVs{length(EVs)+1} = fullfile(outStimuli{ss},stimuli_dirs{j},WAFile{1});
            end
            % Make the EVs for regression
            tmpEV = zeros(TR*size(allTCs,2),1); % time in msec
            allEVs = repmat(tmpEV,1,length(EVs));
            for ee = 1:length(EVs)
                evTimes = 1000*load(EVs{ee}); % convert to msec
                for bb = 1:size(evTimes,1)
                    allEVs(evTimes(bb,1):evTimes(bb,1)+evTimes(bb,2),ee) = 1;
                end
            end
            downEVs = downsample(allEVs,TR);
            downEVs = downEVs - repmat(mean(downEVs),size(downEVs,1),1); % mean center
            outEVs = [ones(size(allTCs,2),1),downEVs];
            % Run the regression
            betaWeights = outEVs\allTCs';
        end
    end
end