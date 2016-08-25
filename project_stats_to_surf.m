function project_stats_to_surf(session_dir,subject_name,sentDirName,locDirName,sentNum,locNum,sentOutDir,locOutDir,sentNames,locNames)

%% Projects sentence and localizer volumetric data to the cortical surface
%
%   Usage:
%   project_stats_to_surf(session_dir,subject_name,sentDirName,locDirName,sentNum,locNum,sentOutDir,locOutDir)
%session_dir = '/Users/Shared/ClusterCopy/Data/AMMet/102MM/M091515M'
% subject_name = 'M091515M'
% sentDirName = 'Series_003_BOLD_METFAM1'
% locDirName = 'Series_011_BOLD_LOCAUD' 
% sentNum = 'sent1aud.feat'
% locNum = 'locaud.feat'
% sentOutDir = '/Users/Shared/ClusterCopy/Data/AMMet/102MM/M091515M/test1'
% locOutDir = '/Users/Shared/ClusterCopy/Data/AMMet/102MM/M091515M/test2'
%   tksurfer L092915B lh inflated -overlay /Users/Shared/ClusterCopy/Data/AMMet/104LB/L092915B/locout/lh.zstat4.nii.gz -overlay /Users/Shared/ClusterCopy/Data/AMMet/104LB/L092915B/sentout/lh.zstat4.nii.gz
%

%% set defaults
% stat files
sentDir = fullfile(session_dir,sentDirName);
locdir = fullfile(session_dir,locDirName);
sentStatDir = fullfile(session_dir,sentNum,'stats');
sentMeanDir = fullfile(session_dir,sentNum);
locStatDir = fullfile(session_dir,locNum,'stats');
% reg files
sentregFile = fullfile(sentDir,'func_bbreg.dat');
locregFile = fullfile(locdir,'func_bbreg.dat');
% out names
hemis = {'lh' 'rh'};
if ~exist('sentNames','var')
    sentNames = {...
        %'zstat1'...
        %'zstat2'...
        %'zstat3' ...
        %'zstat4' ...
        %'cope1'...
        %'cope2'...
        %'cope3'...
        %'cope4'...
        'mean_func'...
        };
end
if ~exist('locNames','var')
    locNames = {...
        %'zstat1' ...
        %'zstat2' ...
        %'zstat3' ...
        %'zstat4' ...
        %'zstat5' ...
        %'zstat6' ...
        %'zstat7' ...
        %'zstat8' ...
        %'zstat9' ...
        %'zfstat1'...
        %s'zfstat2'...
        };
end
if ~exist(sentOutDir,'dir')
    mkdir(sentOutDir);
end
if ~exist(locOutDir,'dir')
    mkdir(locOutDir);
end
%% Project sentences
for hh = 1:length(hemis)
    for oo = 1:length(sentNames)
        hemi = hemis{hh};
        if ~strcmp(sentNames{oo},'mean_func')
            tmpin = fullfile(sentStatDir,[sentNames{oo} '.nii.gz']);
            tmpout = fullfile(sentOutDir,[hemi '.' sentNames{oo} '.nii.gz']);
            system(['mri_vol2surf --mov ' tmpin ' --reg ' sentregFile ...
                ' --hemi ' hemi ' --projfrac-avg 0 1 0.1 --o ' tmpout]);
        else
            tmpin = fullfile(sentMeanDir,[sentNames{oo} '.nii.gz']);
            tmpout = fullfile(sentOutDir,[hemi '.' sentNames{oo} '.nii.gz']);
            system(['mri_vol2surf --mov ' tmpin ' --reg ' sentregFile ...
                ' --hemi ' hemi ' --projfrac-avg 0 1 0.1 --o ' tmpout]);
        end
    end
end
%% Project localizers
for hh = 1:length(hemis)
    for oo = 1:length(locNames)
        hemi = hemis{hh};
        tmpin = fullfile(locStatDir,[locNames{oo} '.nii.gz']);
        tmpout = fullfile(locOutDir,[hemi '.' locNames{oo} '.nii.gz']);
        system(['mri_vol2surf --mov ' tmpin ' --reg ' locregFile ...
            ' --hemi ' hemi ' --projfrac-avg 0 1 0.1 --o ' tmpout]);
    end
end
% %% view the surface output
% for hh = 1:length(hemis)
%     for oo = 1:length(sentNames)
%         hemi = hemis{hh};
%         tmpind = fullfile(sentOutDir,[hemi '.' sentNames{oo} '.nii.gz']);
%         thresh = load_nifti(tmpind); % load the surface file
%         surface_plot('zstat',tmpind,subject_name,hemi,'inflated',thresh.vol>2)
%     end
% end
% %% view the surface output
% for hh = 1:length(hemis)
%     for oo = 1:length(locNames)
%         hemi = hemis{hh};
%         tmpind = fullfile(locOutDir,[hemi '.' locNames{oo} '.nii.gz']);
%         thresh = load_nifti(tmpind); % load the surface file
%         surface_plot('zstat',tmpind,subject_name,hemi,'inflated',thresh.vol>2)
%     end
end