function project_feat2anat(session_dir,subject_name,func,runs,Comps,SUBJECTS_DIR)

% Projects feat stats files to anatomical space, in feat/stats dir
%
%   Usage:
%   project_feat_to_anat(session_dir,subject_name,func,runs,Comps,hemis,SUBJECTS_DIR)
%
%   numComps = number of comparisons (default = 6)
%
%   Written by Andrew S Bock Sep 2015

%% set defaults
if ~exist('Comps','var')
    Comps = 1:6; %  comparisons
end
if ~exist('func','var')
    func = 's5.rf.tf';
end
% get bold runs
b = find_bold(session_dir);
if ~exist('runs','var')
    runs = 1:length(b);
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
stats = {'cope' 'tstat' 'zstat'};
%% Project stats to anatomical space
for i = runs
    tmpreg = listdir(fullfile(session_dir,b{i},'*bbreg.dat'),'files');
    bbreg = fullfile(session_dir,b{i},tmpreg{1});
    for j = Comps
        for z = 1:length(stats)
            outfile = fullfile(session_dir,b{i},[func '.feat'],'stats',...
                [stats{z} num2str(j) '.anat.nii.gz']);
            [~,~] = system(['mri_vol2vol --mov ' ...
                fullfile(session_dir,b{i},[func '.feat'],'stats',...
                [stats{z} num2str(j) '.nii.gz']) ' --targ ' ...
                fullfile(SUBJECTS_DIR,subject_name,'mri','T1.mgz') ' --o ' ...
                outfile ' --reg ' bbreg]);
        end
    end
end