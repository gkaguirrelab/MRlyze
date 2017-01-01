function combine_template_hemispheres(session_dir,subject_name,pRF)

% Projects the 'areas' 'ecc' and 'pol' surface templates to the respective
% bold run volumes, and combines across hemispheres.
%
%   Usage:
%   combine_template_hemispheres(session_dir,subject_name,pRF)
%
%   default:
%   pRF = 1; % assumes 'areas_pRF', 'ecc_pRF', and 'pol_pRF' templates
%   func = 'brf';
%
%   Written by Andrew S Bock Oct 2015

%% set defaults
if ~exist('pRF','var')
    pRF = 1;
end
if pRF
    surfvals = {'areas_pRF' 'ecc_pRF' 'pol_pRF'};
else
    surfvals = {'areas' 'ecc' 'pol'};
end
hemis = {'lh' 'rh'};
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
%% Project surface templates to bold run volumes
for hh = 1:length(hemis)
    hemi = hemis{hh};
    for ss = 1:length(surfvals)
        surfval = fullfile(session_dir,[hemi '.' surfvals{ss} '.nii.gz']);
        tempVol = fullfile(SUBJECTS_DIR,subject_name,'mri/T1.mgz');
        outFile = fullfile(session_dir,[hemi '.' surfvals{ss} '.vol.nii.gz']);
        system(['mri_surf2vol --surfval ' surfval ' --hemi ' hemi ' --fillribbon ' ...
            ' --identity ' subject_name ' --template ' tempVol ' --o ' outFile]);
    end
end
%% Combine hemis
for ss = 1:length(surfvals)
    lhvol = fullfile(session_dir,['lh.' surfvals{ss} '.vol.nii.gz']);
    rhvol = fullfile(session_dir,['rh.' surfvals{ss} '.vol.nii.gz']);
    lh = load_nifti(lhvol);
    rh = load_nifti(rhvol);
    % set nans to zero
    lh.vol(isnan(lh.vol)) = 0;
    rh.vol(isnan(rh.vol)) = 0;
    % set voxels in both hemispheres to zero
    lhind = lh.vol~=0;
    rhind = rh.vol~=0;
    mhind = lhind & rhind; % shared voxels across hemis are excluded
    lh.vol(mhind) = 0;
    rh.vol(mhind) = 0;
    % Combine hemisphere volumes and save
    mh = lh;
    mh.vol = lh.vol + rh.vol;
    save_nifti(mh,fullfile(session_dir,['mh.' surfvals{ss} '.vol.nii.gz']));
end