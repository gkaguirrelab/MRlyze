function template_pRF(session_dir,subject_name,pRF_dir,data_file)

% Takes the output from Mathematica template fitting (Benson, 2012;2014),
%   seperates into individual maps, transforms back to .nii.gz, and
%   converts polar angles back to radians.
%
% Written by Andrew S Bock Feb 2015

%% Set defaults
if ~exist('pRF_dir','var')
    pRF_dir = fullfile(session_dir,'prfs');
end
if ~exist('data_file','var')
    data_file = fullfile(pRF_dir,'template_pRFs.mgz');
end

%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,subject_name,pRF_dir,data_file)

%% Load in data
disp('Loading data...');
mgz = load_mgh(data_file);mgz = squeeze(mgz);
mgz(mgz==99999) = nan;
% Load in a nifti file to overwrite with template values
nii = load_nifti(fullfile(pRF_dir,'mh.ecc.avg_sym.nii.gz'));
disp('done.');
%% Save eccentricity, polar angle, and correlation values
% note these correlation values relate to the pRFs, not to the template fit
disp('Saving maps as nifti files...');
nii.vol = mgz(:,1);
nii.vol = deg2rad((nii.vol - 90)); % convert deg to radians
save_nifti(nii,fullfile(pRF_dir,'mh.pol.avg_template.nii.gz'));
nii.vol = mgz(:,2);
save_nifti(nii,fullfile(pRF_dir,'mh.ecc.avg_template.nii.gz'));
nii.vol = mgz(:,3);
save_nifti(nii,fullfile(pRF_dir,'mh.co.avg_template.nii.gz'));
disp('done.');
%% Project to the subject's native fit to fsaverage_sym space
disp('Projecting to fsaverage_sym space...');
maps = {'ecc' 'pol' 'co'};
for m = 1:length(maps)
    savename = fullfile(pRF_dir,['mh.' maps{m} '.avg']);
    [~,~] = system(['mri_surf2surf --hemi lh --srcsubject ' subject_name ...
        ' --srcsurfval ' savename '_template.nii.gz --trgsubject fsaverage_sym' ...
        ' --trgsurfval ' savename '_sym_template.nii.gz']);
end
disp('done.');
