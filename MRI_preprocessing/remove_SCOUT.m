function remove_SCOUT(session_dir)

%   Removes the SCOUT directory, as single TR BOLD run typically run at the
%   beginning of a MRI session. This function removes any directory in the
%   session_dir whose name contains the string 'SCOUT'.
%
%   Usage: remove_SCOUT(session_dir)
%
%   Written by Andrew S Bock Apr 2015
%

%% Add to log
SaveLogInfo(session_dir, mfilename,session_dir)
%% Remove SCOUT directory
d=listdir(session_dir,'dirs');
scoutDir = d(~cellfun(@isempty,strfind(d,'SCOUT')));
if ~isempty(scoutDir);
    for s = 1:length(scoutDir);
        disp(['Removing SCOUT directory: ' fullfile(session_dir,scoutDir{s})]);
        rmdir(fullfile(session_dir,scoutDir{s}), 's');
    end
else
    disp('No SCOUT directory found.');
end
disp('done.');