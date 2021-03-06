function [val] = get_HRF_peak(session_dir,hemi,srcROI,trgROI,map_type,ROI)

% Returns the mean HRF peak for a given ROI
%
%   Usage:
%   [val] = get_HRF_peak(session_dir,srcROI,trgROI,map_type,ROI)
%
%   Written by Andrew S Bock Jul 2015

%%
prfpol = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.copol.prfs.nii.gz']));
prfsig = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.cosig1.prfs.nii.gz']));
prfco = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.co.prfs.nii.gz']));
prfpeak = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.copeakt.prfs.nii.gz']));

cfpeak = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.cotrgpeakt.cfs.nii.gz']));
cfshift = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.coshiftt.cfs.nii.gz']));

if strcmp(srcROI,'volume')
    mask = load_nifti(fullfile(session_dir,[hemi '.' ROI '.mask.nii.gz']));
    maskind = (mask.vol > 0);
    val.prfmask = mean(prfpeak.vol(maskind));
    val.cfmask = mean(cfpeak.vol(maskind));
    val.cfmask_shift = mean(cfshift.vol(maskind));
elseif strcmp(srcROI,'cortex');
    if strcmp(hemi,'lh') || strcmp(hemi,'rh')
        mask = load_nifti(fullfile(session_dir,[hemi '.areas_pRF.nii.gz']));
    elseif strcmp(hemi,'mh')
        mask = load_nifti(fullfile(session_dir,'lh.areas_pRF.nii.gz'));        
    end
    if strcmp(ROI,'V1')
        maskind = ((mask.vol == -1 | mask.vol == 1) & prfco.vol > 0.3873 & prfsig.vol > 0.25);
        dind = (mask.vol == -1 & prfco.vol > 0.3873 & prfsig.vol > 0.25);
        vind = (mask.vol == 1 & prfco.vol > 0.3873 & prfsig.vol > 0.25);
    elseif strcmp(ROI,'V2')
        maskind = ((mask.vol == -2 | mask.vol == 2) & prfco.vol > 0.3873 & prfsig.vol > 0.25);
        dind = (mask.vol == -2 & prfco.vol > 0.3873 & prfsig.vol > 0.25);
        vind = (mask.vol == 2 & prfco.vol > 0.3873 & prfsig.vol > 0.25);
        if strcmp(hemi,'lh') || strcmp(hemi,'mh')
            dpolind = prfpol.vol > 0 & prfpol.vol < pi/2;
            vpolind = prfpol.vol < 0 & prfpol.vol > -pi/2;
        elseif strcmp(hemi,'rh')
            dpolind = prfpol.vol > pi/2 & prfpol.vol < pi;
            vpolind = prfpol.vol < -pi/2 & prfpol.vol > -pi;
        end
        dind(~dpolind) = 0;
        vind(~vpolind) = 0;
    elseif strcmp(ROI,'V3')
        maskind = ((mask.vol == -3 | mask.vol == 3) & prfco.vol > 0.3873 & prfsig.vol > 0.25);
        dind = (mask.vol == -3 & prfco.vol > 0.3873 & prfsig.vol > 0.25);
        vind = (mask.vol == 3 & prfco.vol > 0.3873 & prfsig.vol > 0.25);
        if strcmp(hemi,'lh') || strcmp(hemi,'mh')
            dpolind = prfpol.vol > 0 & prfpol.vol < pi/2;
            vpolind = prfpol.vol < 0 & prfpol.vol > -pi/2;
        elseif strcmp(hemi,'rh')
            dpolind = prfpol.vol > pi/2 & prfpol.vol < pi;
            vpolind = prfpol.vol < -pi/2 & prfpol.vol > -pi;
        end
        dind(~dpolind) = 0;
        vind(~vpolind) = 0;
    end
    val.prfmask = mean(prfpeak.vol(maskind));
    val.prfd = mean(prfpeak.vol(dind));
    val.prfv = mean(prfpeak.vol(vind));
    val.cfmask = mean(cfpeak.vol(maskind));
    val.cfd = mean(cfpeak.vol(dind));
    val.cfv = mean(cfpeak.vol(vind));
    val.cfmask_shift = mean(cfshift.vol(maskind));
    val.cfd_shift = mean(cfshift.vol(dind));
    val.cfv_shift = mean(cfshift.vol(vind));
end

