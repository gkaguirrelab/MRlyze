function project_MP_ROIs_2func(session_dir,runNums,func)

% Projects M and P ROIs to functional space
%
%   Usage:
%   project_MP_ROIs_2func(session_dir,runNums)
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
hemis = {'lh' 'rh'};
ROIs = {'LGN' 'V1_pRF' 'V2_pRF' 'V3_pRF'};
%% Find bold dirs
b = find_bold(session_dir);

%% Create matrices for each run
for hh = 1:length(hemis)
    hemi = hemis{hh};
    for rr = runNums
        for ro = 1:length(ROIs)
            ROI = ROIs{ro};
            % Project ROI to functional space
            targvol = fullfile(session_dir,b{rr},'drf.tf.nii.gz');
            reg = fullfile(session_dir,b{rr},'rf_bbreg.dat');
            Minvol = fullfile(session_dir,[hemi '.' func '.' ROI '.M.nii.gz']);
            Moutvol = fullfile(session_dir,b{rr},[hemi '.' func '.' ROI '.M.nii.gz']);
            Pinvol = fullfile(session_dir,[hemi '.' func '.' ROI '.P.nii.gz']);
            Poutvol = fullfile(session_dir,b{rr},[hemi '.' func '.' ROI '.P.nii.gz']);
            mri_vol2vol(targvol,Minvol,Moutvol,reg,'nearest',1);
            mri_vol2vol(targvol,Pinvol,Poutvol,reg,'nearest',1);
        end
    end
end