function rerun_regress(inScript,scriptDir,logDir)

% reruns regress_template scripts on the cluster. Assume the log files in
% the 'logDir' are ONLY from 'inScript'
%
%   Usage:
%   rerun_regress(inScript,scriptDir,logDir)
%
%   Written by Andrew S Bock Mar 2016

%% Set defaults
if ~exist('LOG_dir','var')
    logDir = '/data/jet/abock/LOGS';
end
thisDir = pwd;
%% Load inputScript
ct = 0;
fid = fopen(fullfile(scriptDir,inScript));
while 1
    ct = ct + 1;
    tmp = fgetl(fid);
    if ~ischar(tmp), break, end
    shScripts{ct} = tmp;
end
fclose(fid);
%% Get error files
e = listdir(fullfile(logDir,'./*.e*'),'files');
for i = 1:length(e)
    fid = fopen(fullfile(logDir,e{i}));
    tmp = fread(fid);
    if ~isempty(tmp)
        cd(scriptDir);
        system(shScripts{i});
    end
    fclose(fid);
end
cd(thisDir);