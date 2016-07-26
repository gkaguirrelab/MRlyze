function clean_up(session_dir)

% Cleans up temporary files following the MRI preprocessing scripts.
%   
% WARNING: only the files needed for the following processing steps are kept,
% everything else is deleted.
%
% Written by Andrew S Bock Jan 2015
% Updated on July 2016 by ASB, GF 

%% define variable names to keep in each bold dir
filesToKeep = { ...
    's5.wdrf.tf.nii.gz' ...
    's5.wdrf.tf.surf.rh.nii.gz' ...
    's5.wdrf.tf.surf.lh.nii.gz' ...
    'wdrf.tf.nii.gz' ...
    };
foldersToKeep = { ...
    'wdrf.tf.feat' ...
    's5.wdrf.tf.feat' ...
    };

%% find bold directories
b = find_bold(session_dir);
%% remove all files and folders other than those to keep
for nn = 1:length(b)
    % check all files
    allFiles = listdir(fullfile(session_dir,b{nn}),'files');
    toKeep = ismember(allFiles,filesToKeep);
    for kk = 1:length(toKeep)
        if toKeep(kk) == 0
            fprintf ('\n Removing %s', fullfile(session_dir,b{nn},allFiles{kk}))
            system(['rm -rf ' fullfile(session_dir,b{nn},allFiles{kk})])
        end
    end
    clear ('toKeep');
    % check all folders
    allFolders = listdir(fullfile(session_dir,b{nn}),'dirs');
    toKeep = ismember(allFolders,foldersToKeep);
    for jj = 1:length(toKeep)
        if toKeep(jj) == 0
            fprintf ('\n Removing %s', fullfile(session_dir,b{nn},allFolders{jj}))
            system(['rm -rf ' fullfile(session_dir,b{nn},allFolders{jj})])
        end
    end
    clear ('toKeep');
end

