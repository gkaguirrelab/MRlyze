function [F,df] = run_F_test(session_dir,subject_name,func,runs,Comps,hemis,SUBJECTS_DIR)

% Runs 'F_test' on feat data
%
%   Usage:
%   [F,df] = run_F_test(session_dir,subject_name,func,runs,Comps,hemis,SUBJECTS_DIR)
%
%   Output (in session_directory):
%       <func>.Fstat.anat.nii.gz - Fstats in anatomical volume
%       <func>.Fstat.surf.lh.nii.gz - Fstats on surface (left hemisphere)
%       <func>.Fstat.surf.rh.nii.gz - Fstats on surface (right hemisphere)
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
if ~exist('hemis','var')
    hemis = {'lh' 'rh'};
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
%% Project copes to anatomical space
for i = runs
    tmpreg = listdir(fullfile(session_dir,b{i},'*bbreg.dat'),'files');
    bbreg = fullfile(session_dir,b{i},tmpreg{1});
    for j = Comps
        outfile = fullfile(session_dir,b{i},[func '.feat'],'stats',...
            ['cope' num2str(j) '.anat.nii.gz']);
        [~,~] = system(['mri_vol2vol --mov ' ...
            fullfile(session_dir,b{i},[func '.feat'],'stats',...
            ['cope' num2str(j) '.nii.gz']) ' --targ ' ...
            fullfile(SUBJECTS_DIR,subject_name,'mri','T1.mgz') ' --o ' ...
            outfile ' --reg ' bbreg]);
    end
end
%% get copes
mat = nan(length(runs),length(Comps),256,256,256);
ct = 0;
for i = runs
    ct = ct + 1;
    jct = 0;
    for j = Comps
        jct = jct + 1;
        tmp = load_nifti(fullfile(session_dir,b{i},[func '.feat'],'stats',...
            ['cope' num2str(j) '.anat.nii.gz']));
        mat(ct,jct,:,:,:) = tmp.vol;
    end
end
%% Run F test
matdims = size(mat);
flatmat = reshape(mat,matdims(1),matdims(2),matdims(3)*matdims(4)*matdims(5));
[F,df] = F_test(flatmat);
%% Save output
F = reshape(F,matdims(3),matdims(4),matdims(5));
tmp.vol = F;
Ffile = fullfile(session_dir,[func '.Fstat.anat.nii.gz']);
save_nifti(tmp,Ffile);
for hh = 1:length(hemis)
    hemi = hemis{hh};
    outfile = fullfile(session_dir,[func '.Fstat.surf.' hemi '.nii.gz']);
    [~,~] = system(['mri_vol2surf --mov ' ...
        Ffile ' --regheader ' subject_name ' --hemi ' hemi ' --o ' outfile]);
end
