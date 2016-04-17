function clean_up(session_dir)

% Cleans up temporary files following the MRI preprocessing scripts.
%   Files/directories removed in the session_dir include:
%
%   <bold_dir>/*.feat
%   <bold_dir>/brf.*
%   <bold_dir>/brf_spikes*
%   <bold_dir>/f.nii.gz
%   <bold_dir>/raw_f.nii.gz
%   aseg.*
%   timestamp
%   parfor_progress
%
% Note - <bold_dir> could be /*BOLD_*, '/*bold*, /*EPI_*, /*ep2d* or 'RUN*'
%
% Written by Andrew S Bock Jan 2015

%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir);

%% Find bold run directories
d = listdir(fullfile(session_dir,'*BOLD_*'),'dirs');
if isempty(d)
    d = listdir(fullfile(session_dir,'*bold_*'),'dirs');
end
if isempty(d)
    d = listdir(fullfile(session_dir,'*EPI_*'),'dirs');
end
if isempty(d)
    d = listdir(fullfile(session_dir,'*ep2d*'),'dirs');
end
if isempty(d)
    d = listdir(fullfile(session_dir,'RUN*'),'dirs');
end
nruns = length(d);
disp(['Session_dir = ' session_dir]);
disp(['Number of runs = ' num2str(nruns)]);
%% Remove files and directories
for r = 1:nruns
    % Remove feat directories
    featdirs = listdir(fullfile(session_dir,d{r}),'dirs');
    if ~isempty(featdirs)
        for f = 1:length(featdirs)
            if strcmp(featdirs{f}(end-4:end),'.feat')
                system(['rm -r ' fullfile(session_dir,d{r},featdirs{f})]);
            end
        end
    end
    % Remove intermediate files in bold directories
    brf = listdir(fullfile(session_dir,d{r},'brf.*'),'files');
    if ~isempty(brf)
        for b = 1:length(brf)
            delete(fullfile(session_dir,d{r},brf{b}))
        end
    end
    spikes = listdir(fullfile(session_dir,d{r},'brf_spikes*'),'files');
    if ~isempty(spikes)
        for s = 1:length(spikes)
            delete(fullfile(session_dir,d{r},spikes{s}));
        end
    end
    if exist(fullfile(session_dir,d{r},'f.nii.gz'),'file')
        delete(fullfile(session_dir,d{r},'f.nii.gz'));
    end
    if exist(fullfile(session_dir,d{r},'raw_f.nii.gz'),'file')
        delete(fullfile(session_dir,d{r},'raw_f.nii.gz'));
    end
end
% Remove intermediate files in session directory
aseg = listdir(fullfile(session_dir,'aseg.*'),'files');
if ~isempty(aseg)
    for a = 1:length(aseg)
        delete(fullfile(session_dir,aseg{a}));
    end
end