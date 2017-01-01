function make_covariates_7T_OneLight(session_dir,runNum,mat_num,mat_dir)

% Outputs covariate text files for FSL's FEAT
%
%   Usage:
%   make_covariates_7T_OneLight(session_dir,runNum,mat_num,mat_dir)
%
%   session_dir = session directory
%   runNum = bold directory in session_directory
%   mat_num = mat file corresponding to bold directory
%   mat_dir = directory containing mat files from stimulus computer
%       default = fullfile(session_dir,'OneLightMatFiles');
%
%   NOTE: assumes the first flicker frequency is 0Hz, so the condition text
%   files that are written start with the 2nd flicker index value
%
%   Written by Andrew S Bock Aug 2015

%% Set defaults
if ~exist('matdir','var')
    mat_dir = fullfile(session_dir,'OneLightMatFiles');
end
d = find_bold(session_dir);
m = listdir(mat_dir,'files');
out_dir = fullfile(session_dir,d{runNum});
mat_file = fullfile(mat_dir,m{mat_num});
data = load(mat_file);
%% Get basic trial information
numTrials = data.params.nTrials; % number of blocks
startRun = data.params.responseStruct.tBlockStart;
endRun = data.params.responseStruct.tBlockEnd;
lengthRun = endRun - startRun;
%% Find attention tasks
attTask = zeros(numTrials,1);
for i = 1:numTrials
    % find if there was an attention task
    if data.params.responseStruct.events(i).attentionTask.segmentFlag
        attTask(i) = 1;
    end
end
%% Find attention task trials where subject didn't push a button
ct = 0;
for i = 1:numTrials
    % find if there was an attention task
    if attTask(i)
        % find if the subject did NOT push a button during attention task
        if isempty(data.params.responseStruct.events(i).buffer)
            ct = ct + 1;
            startTrial = data.params.responseStruct.events(i).tTrialStart;
            endTrial = data.params.responseStruct.events(i).tTrialEnd;
            badTrials(ct,1) = startTrial - startRun;
            badTrials(ct,2) = endTrial - startTrial;
            badTrials(ct,3) = 1;
        end
    end
end
%% Find timings of attention task flashes
taskTrials = [0 0 0];
ct = 0;
for i = 1:numTrials
    % find if there was an attention task
    if attTask(i)
        ct = ct + 1;
        % start time (1) (in seconds from start of run)
        startTask = data.params.responseStruct.events(i).t(...
            data.params.responseStruct.events(i).attentionTask.T == 1) - ...
            data.params.responseStruct.tBlockStart;
        taskTrials(ct,1) = startTask;
        % duration (in seconds from start of run)
        endTask = data.params.responseStruct.events(i).t(...
            data.params.responseStruct.events(i).attentionTask.T == -1) - ...
            data.params.responseStruct.tBlockStart;
        taskTrials(ct,2) = endTask - startTask;
        % fill in third column with ones ("strength" of the event)
        taskTrials(ct,3) = 1;
    end
end
%% Find timings of various frequencies
FreqOrder = data.params.theFrequencyIndices; % Frequency indices
Freqs = unique(FreqOrder); % Get unique frequencies
for i = 1:length(Freqs)
    ct = 0;
    for j = 1:length(FreqOrder)
        if Freqs(i) == FreqOrder(j)
            ct = ct + 1;
            startTrial = data.params.responseStruct.events(j).tTrialStart;
            endTrial = data.params.responseStruct.events(j).tTrialEnd;
            contrasts(ct,1,i) = startTrial - startRun;
            contrasts(ct,2,i) = endTrial - startTrial;
            contrasts(ct,3,i) = 1;
        end
    end
end
%% Write out the covariate files
if exist('badTrials','var')
    dlmwrite(fullfile(out_dir,'bad_trials.txt'),badTrials,'delimiter',' ','precision','%10.5f');
end
dlmwrite(fullfile(out_dir,'task_trials.txt'),taskTrials,'delimiter',' ','precision','%10.5f');
for i = 1:size(contrasts,3)
    if i>1 % first condition is 0Hz flicker
        dlmwrite(fullfile(out_dir,['condition_' num2str(i-1) '.txt']),contrasts(:,:,i),'delimiter',' ','precision','%10.5f')
    end
end