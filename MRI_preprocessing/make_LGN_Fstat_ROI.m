function make_LGN_Fstat_ROI(session_dir,subject_name,runNums,func,thresh)

% Makes a mask in the LGN based on anatomical ROI and Fstat values
%
%   Usage:
%
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
if ~exist('thresh','var')
    thresh = 0.05; % p-value threshold
end
hemis = {'lh' 'rh'};
%% get bold dirs
d = find_bold(session_dir);

%% Load Fstat values from runs
for i = runNums
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        convert_F_to_p(session_dir,subject_name,runNums,func,thresh);
        statsDir = fullfile(session_dir,d{i},[func '.feat'],'stats');
        pmask = load_nifti(fullfile(statsDir,'pval.mask.anat.nii.gz'));
        LGNmask = load_nifti(fullfile(session_dir,[hemi '.LGN.nii.gz']));
        dotmask = load_nifti(fullfile(session_dir,[hemi '.dots.LGN.nii.gz']));
    end
end