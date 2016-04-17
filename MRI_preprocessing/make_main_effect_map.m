function make_main_effect_map(session_dir,subject_name,func,contrastNum,hemis,SUBJECTS_DIR)

% Makes volume and surface images representing the main effect, as found in
% a feat directory (specified by 'contrast_Num' input)
%
%   Usage:
%   make_main_effect_map(session_dir,subject_name,featdir,contrastNum,hemis,SUBJECTS_DIR)
%
%   Output (in session_directory):
%       <featdir>.main_effect.anat.nii.gz - main effect in anatomical volume
%       <featdir>.main_effect.surf.lh.nii.gz - main effect on surface (left hemisphere)
%       <featdir>.main_effect.surf.rh.nii.gz - main effect on surface (right hemisphere)
%
%   Written by Andrew S Bock Sep 2015

%% set defaults
if ~exist('func','var')
    func = 'dbrf.tf';
end
if ~exist('hemis','var')
    hemis = {'lh' 'rh'};
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
%% get bold runs
b = find_bold(session_dir);
%% Project zstats to anatomical space
for i = 1:length(b)
    outfile = fullfile(session_dir,b{i},[func '.feat'],'stats',...
        ['zstat' num2str(contrastNum) '.anat.nii.gz']);
    [~,~] = system(['mri_vol2vol --mov ' ...
        fullfile(session_dir,b{i},[func '.feat'],'stats',...
        ['zstat' num2str(contrastNum) '.nii.gz']) ' --targ ' ...
        fullfile(SUBJECTS_DIR,subject_name,'mri','T1.mgz') ' --o ' ...
        outfile ' --reg ' ...
        fullfile(session_dir,b{i},'brf_bbreg.dat')]);
end
%% get zstats
mat = nan(length(b),256,256,256);
for i = 1:length(b)
    tmp = load_nifti(fullfile(session_dir,b{i},[func '.feat'],'stats',...
        ['zstat' num2str(contrastNum) '.anat.nii.gz']));
    mat(i,:,:,:) = tmp.vol;
end
%% Average runs to get mean main effect
main_effect = squeeze(mean(mat));
%% Save output
tmp.vol = main_effect;
mfile = fullfile(session_dir,[func '.main_effect.anat.nii.gz']);
save_nifti(tmp,mfile);
for hh = 1:length(hemis)
    hemi = hemis{hh};
    outfile = fullfile(session_dir,[func '.main_effect.surf.' hemi '.nii.gz']);
    [~,~] = system(['mri_vol2surf --mov ' ...
        mfile ' --regheader ' subject_name ' --hemi ' hemi ' --o ' outfile]);
end
