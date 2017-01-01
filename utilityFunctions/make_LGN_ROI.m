function make_LGN_ROI(session_dir,subject_name,thresh,SUBJECTS_DIR)

% An LGN ROI was created in Freesurfer's cvs_MNI space, based on FSL's 1mm
% atlas. This function then takes that ROI and projects it to subject
% native space.
%
%   Usage:
%   make_LGN_ROI(session_dir,subject_name,thresh,SUBJECTS_DIR)
%
%   Outputs (in 'session_dir/anat_templates'):
%       lh.LGN.prob.nii.gz - voxel values reflect the probability of being
%           in the left LGN
%       rh.LGN.prob.nii.gz - voxel values reflect the probability of being
%           in the right LGN
%       lh.LGN.nii.gz - binary mask volume, threshold using the input
%           'thresh' <default = 25>
%       rh.LGN.nii.gz - binary mask volume, threshold using the input
%           'thresh' <default = 25>
%
%   Written by Andrew S Bock Sep 2015

%% set defaults
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
if ~exist('thresh','var')
    thresh = 25; % probability of LGN in native space
end
outDir = fullfile(session_dir,'anat_templates');
if ~exist(outDir,'dir')
    mkdir(outDir);
end
in_vols = {...
    '~/data/lh.LGN.prob.nii.gz' ...
    '~/data/rh.LGN.prob.nii.gz' ...
    };
prob_vols = {...
    fullfile(outDir,'lh.LGN.prob.nii.gz') ...
    fullfile(outDir,'rh.LGN.prob.nii.gz') ...
    };
out_vols = {...
    fullfile(outDir,'lh.LGN.nii.gz') ...
    fullfile(outDir,'rh.LGN.nii.gz') ...
    };
%% Create LGN ROIs for left and right hemisphere
for i = 1:length(in_vols)
    in_vol = in_vols{i};
    prob_vol = prob_vols{i};
    ref_vol = fullfile(SUBJECTS_DIR,subject_name,'mri','T1.mgz');
    apply_cvs_inverse(in_vol,prob_vol,ref_vol,subject_name,SUBJECTS_DIR)
end
%% Threshold voxels
for i = 1:length(out_vols)
    prob_vol = prob_vols{i};
    out_vol = out_vols{i};
    tmp = load_nifti(prob_vol);
    tmp.vol(tmp.vol<thresh) = 0;
    tmp.vol(tmp.vol>0) = 1;
    save_nifti(tmp,out_vol);
end
%% Combine across hemispheres
lh = load_nifti(out_vols{1});
rh = load_nifti(out_vols{2});
mh = lh;
mh.vol = lh.vol + rh.vol;
mh.vol(mh.vol>0) = 1;
save_nifti(mh,fullfile(outDir,'mh.LGN.nii.gz'));
