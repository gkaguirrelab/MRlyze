function create_MP_correlation_matrix(session_dir,runNums,tcfunc,roifunc,ROIs,hemis)

% Creates a cross-correlation matrix between specified 'ROIs' in each
% hemissphere, for each run specified by 'runNums'.
%
%   Usage:
%
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
if ~exist('ROIs','var')
    ROIs = {'LGN.M' 'V1_pRF.M' 'V2_pRF.M' 'V3_pRF.M' ...
        'LGN.P' 'V1_pRF.P' 'V2_pRF.P' 'V3_pRF.P' ...
        };
end
if ~exist('hemis','var')
    hemis = {'lh' 'rh'};
end
%% find bold directories
b = find_bold(session_dir);
%% Make correlation matrices
for hh = 1:length(hemis)
    hemi = hemis{hh};
    for rr = runNums
        boldDir = fullfile(session_dir,b{rr});
        tc = load_nifti(fullfile(boldDir,[tcfunc '.nii.gz']));
        dims = size(tc.vol);
        flattc = (reshape(tc.vol,dims(1)*dims(2)*dims(3),dims(4)))';
        lhrotc = nan(dims(4),length(ROIs));
        for ro = 1:length(ROIs)
            ROI = ROIs{ro};
            lhtmp = load_nifti(fullfile(boldDir,[hemi '.' roifunc '.' ROI '.nii.gz']));
            lhroInd = lhtmp.vol>0;
            lhroInd = reshape(lhroInd,dims(1)*dims(2)*dims(3),1);
            lhrotc(:,ro) = mean(flattc(:,lhroInd),2);
        end
        corrMat = corr(lhrotc);
        save(fullfile(boldDir,[hemi '.TCvol.' tcfunc '.ROIvol.' roifunc '.mat']),'corrMat','ROIs');
    end
end
%% Make correlation matrices between hemisphers
for rr = runNums
    boldDir = fullfile(session_dir,b{rr});
    tc = load_nifti(fullfile(boldDir,[tcfunc '.nii.gz']));
    dims = size(tc.vol);
    flattc = (reshape(tc.vol,dims(1)*dims(2)*dims(3),dims(4)))';
    for ro = 1:length(ROIs)
        ROI = ROIs{ro};
        lhtmp = load_nifti(fullfile(boldDir,['lh.' roifunc '.' ROI '.nii.gz']));
        rhtmp = load_nifti(fullfile(boldDir,['rh.' roifunc '.' ROI '.nii.gz']));
        lhroInd = lhtmp.vol>0;
        rhroInd = rhtmp.vol>0;
        lhroInd = reshape(lhroInd,dims(1)*dims(2)*dims(3),1);
        rhroInd = reshape(rhroInd,dims(1)*dims(2)*dims(3),1);
        lhrotc = mean(flattc(:,lhroInd),2);
        rhrotc = mean(flattc(:,rhroInd),2);
        mhcorrMat(ro) = corr(lhrotc,rhrotc);
    end
    save(fullfile(boldDir,['mh.TCvol.' tcfunc '.ROIvol.' roifunc '.mat']),'mhcorrMat','ROIs');
end