function project_LGN_MT(session_dir,subject_name,thresh,template_files,hemis,SUBJECTS_DIR)

% Projects anatomical template files in cvs MNI space to subject
%   native space, using FSL's flirt and fnirt
%
%   Usage: project_LGN_MT(session_dir,subject_name,template_files,hemis,thresh,SUBJECTS_DIR,FSLDIR)
%
%   Written by Andrew S Bock Feb 2015

%% Set defaults
if ~exist('template_files','var')
    template_files = {...
        '~/data/lh.LGN.prob.nii.gz' ...
        '~/data/rh.LGN.prob.nii.gz' ...
        '~/data/lh.MT.prob.nii.gz' ...
        '~/data/rh.MT.prob.nii.gz' ...
        };
end
if ~exist('thresh','var')
    thresh = 25;
end
if ~exist('hemis','var');
    hemis = {'lh' 'rh'};
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
maps = {'LGN' 'MT'};
%% Setup files
% T1 files
MNIT1 = fullfile(SUBJECTS_DIR,'cvs_avg35_inMNI152/mri/T1.mgz');
SUBJT1 = fullfile(SUBJECTS_DIR,subject_name,'mri/brain.nii.gz');
flirtout = fullfile(SUBJECTS_DIR,subject_name,'mri/MNI_flirt_brain.nii.gz');
fnirtout = fullfile(SUBJECTS_DIR,subject_name,'mri/MNI_fnirt_brain.nii.gz');
flirtmat = fullfile(SUBJECTS_DIR,subject_name,'mri/MNI2brain.mat');
fnirtwarp = fullfile(SUBJECTS_DIR,subject_name,'mri/MNI2brain.nii.gz');
% Output template files
out_vol_files = {...
    fullfile(session_dir,'lh.LGN.prob.nii.gz') ...
    fullfile(session_dir,'rh.LGN.prob.nii.gz') ...
    fullfile(session_dir,'lh.MT.prob.nii.gz') ...
    fullfile(session_dir,'rh.MT.prob.nii.gz') ...
    };
thresh_out_vol_files = {...
    fullfile(session_dir,'lh.LGN.nii.gz') ...
    fullfile(session_dir,'rh.LGN.nii.gz') ...
    fullfile(session_dir,'lh.MT.nii.gz') ...
    fullfile(session_dir,'rh.MT.nii.gz') ...
    };
out_surf_files = {...
    fullfile(session_dir,'lh.MT.prob.surf.nii.gz') ...
    fullfile(session_dir,'rh.MT.prob.surf.nii.gz') ...
    fullfile(session_dir,'lh.MT.surf.nii.gz') ...
    fullfile(session_dir,'rh.MT.surf.nii.gz') ...
    };
%% Register MNI to subject
if ~exist(SUBJT1,'file');
    SUBJmgz = fullfile(SUBJECTS_DIR,subject_name,'mri/brain.mgz');
    system(['mri_convert ' SUBJmgz ' ' SUBJT1]);
end
system(['flirt -in ' MNIT1 ' -ref ' SUBJT1 ' -omat ' flirtmat ' -o ' flirtout]);
%% FNIRT
system(['fnirt --ref=' SUBJT1 ' --in=' MNIT1 ' --aff=' flirtmat ...
    ' --cout=' fnirtwarp ' --iout=' fnirtout]);
%% Project LGN and MT templates to subject volumetric space
disp(['session_dir = ' session_dir]);
disp(['subject_name = ' subject_name]);
disp('Projecting LGN and MT templates to subject space...');
for t = 1:length(template_files)
    system(['applywarp -i ' template_files{t} ' -o ' ...
        out_vol_files{t} ' -r ' SUBJT1 ' -w ' fnirtwarp]);
end
%% Threshold voxels
for i = 1:length(thresh_out_vol_files)
    tmp = load_nifti(out_vol_files{i});
    tmp.vol(tmp.vol<thresh) = 0;
    tmp.vol(tmp.vol>0) = 1;
    save_nifti(tmp,thresh_out_vol_files{i});
end
%% Merge lh and rh into one volume (mh)
for m = 1:length(maps)
    % Probability maps
    lh = load_nifti(fullfile(session_dir,['lh.' maps{m} '.prob.nii.gz']));
    rh = load_nifti(fullfile(session_dir,['rh.' maps{m} '.prob.nii.gz']));
    mhname = fullfile(session_dir,['mh.' maps{m} '.prob.nii.gz']);
    mh = lh;
    mh.vol = zeros(size(mh.vol));
    mh.vol(lh.vol~=0) = lh.vol(lh.vol~=0);
    mh.vol(rh.vol~=0) = rh.vol(rh.vol~=0);
    sharedind = lh.vol~=0 & rh.vol~=0;
    mh.vol(sharedind) = 0; % set voxels shared by both hemis to zero
    save_nifti(mh,mhname);
    % Binary maps
    lh = load_nifti(fullfile(session_dir,['lh.' maps{m} '.nii.gz']));
    rh = load_nifti(fullfile(session_dir,['rh.' maps{m} '.nii.gz']));
    mhname = fullfile(session_dir,['mh.' maps{m} '.nii.gz']);
    mh = lh;
    mh.vol = zeros(size(mh.vol));
    mh.vol(lh.vol~=0) = lh.vol(lh.vol~=0);
    mh.vol(rh.vol~=0) = rh.vol(rh.vol~=0);
    sharedind = lh.vol~=0 & rh.vol~=0;
    mh.vol(sharedind) = 0; % set voxels shared by both hemis to zero
    save_nifti(mh,mhname);
end
%% Project MT to surfae
for hh = 1:length(hemis)
    hemi = hemis{hh};
    system(['mri_vol2surf --hemi ' hemi ' --regheader ' subject_name ' --mov ' ...
        out_vol_files{hh+2} ' --o ' out_surf_files{hh} ' --projfrac-avg 0 1 0.1']);
    system(['mri_vol2surf --hemi ' hemi ' --regheader ' subject_name ' --mov ' ...
        thresh_out_vol_files{hh+2} ' --o ' out_surf_files{hh+2} ' --projfrac 0.5']);
end
disp('done.');