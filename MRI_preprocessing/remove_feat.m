function remove_feat(session_dir,featName)

% deletes feat directories, specified by 'featName' input
%
%   usage:
%   remove_feat(session_dir,featName)
%
%   example:
%   session_dir = '/data/jet/abock/data/Retinotopy/ASB/10272014';
%   featName = 'brf.feat';
%   remove_feat(session_dir,featName)
%
%   written by Anja Jamrozik October 2015
%   updated by Andrew S Bock Oct 2015

%% find bold directories
d = find_bold(session_dir);

%% remove feat directory
for i=1:length(d) 
    disp(['deleting ' fullfile(session_dir,d{i},featName)]);
    system(['rm -rf ' fullfile(session_dir,d{i},featName)]);  
end
disp('done.');