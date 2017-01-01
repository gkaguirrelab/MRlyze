function project_MP_ROIs_2func(session_dir,runNums)

% Creates correlation matrices between ROIs for M and P
%
%   Usage:
%   create_MP_matrices(session_dir,runNums)
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
hemis = {'lh' 'rh'};
ROIs = {'V1_pRF' 'V2_pRF' 'V3_pRF' 'LGN'};
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
            Minvol = fullfile(session_dir,[hemi '.rf.tf.' ROI '.M.nii.gz']);
            Moutvol = fullfile(session_dir,b{rr},[hemi '.rf.tf.' ROI '.M.nii.gz']);
            Pinvol = fullfile(session_dir,[hemi '.rf.tf.' ROI '.P.nii.gz']);
            Poutvol = fullfile(session_dir,b{rr},[hemi '.rf.tf.' ROI '.P.nii.gz']);            
            mri_vol2vol(targvol,Minvol,Moutvol,reg,'nearest',1);
            mri_vol2vol(targvol,Pinvol,Poutvol,reg,'nearest',1);
            M = load_nifti(Moutvol);
            P = load_nifti(Poutvol);
            Mind = find(M.vol>0);
            Pind = find(P.vol>0);
            flicker = load_nifti(fullfile(session_dir,b{rr},'drf.tf.nii.gz'));
            dims = size(flicker.vol);
            flickertc = reshape(flicker.vol,dims(1)*dims(2)*dims(3),dims(4))';
            Mflicker = mean(flickertc(:,Mind));
            Pflicker = mean(flickertc(:,Pind));
        end
    end
end