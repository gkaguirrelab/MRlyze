function project_pRF2vol(session_dir,subject_name)

% Projects the *_pRF.nii.gz retinotopic template maps to the volume
%
%   Usage:
%   project_pRF2vol(session_dir,subject_name)
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
hemis = {'lh' 'rh'};
maps = {'ecc' 'pol' 'areas'};
%% Project maps to volume
for hh = 1:length(hemis)
    hemi = hemis{hh};
    for mm = 1:length(maps)
        map = maps{mm};
        invol = fullfile(session_dir,[hemi '.' map '_pRF.nii.gz']);
        outvol = fullfile(session_dir,[hemi '.' map '_pRF.vol.nii.gz']);
        mri_surf2vol(subject_name,invol,outvol,hemi);
    end
end
