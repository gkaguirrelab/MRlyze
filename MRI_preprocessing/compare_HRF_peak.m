function compare_HRF_peak(session_dir,srcROI,trgROI,map_type,ROIs)

% Creates scatter plots of HRF peadata.
%
%   Usage:
%   compare_CF_pRF(session_dir,srcROI,trgROI,map_type,ROIs)
%
%   Examples:
%   compare_CF_pRF(session_dir,'volume','prf_V1','movie')
%
%   Written by Andrew S Bock Jul 2015
%%
if ~exist('ROIs','var')
    if strcmp(srcROI,'volume');
        ROIs = {'LGN' 'SC'};
    elseif strcmp(srcROI,'cortex');
        ROIs = {'V2' 'V3'};
    end
end
hemis = {'lh' 'rh'};
colors = {'r' 'g' 'b' 'c'};
tboxcoord = {[0.25,0.7,0.1,0.025] [0.25,0.6,0.1,0.025] [0.7,0.7,0.1,0.025] [0.7,0.6,0.1,0.025] ...
    [0.25,0.3,0.1,0.025] [0.25,0.2,0.1,0.025] [0.7,0.3,0.1,0.025] [0.7,0.2,0.1,0.025]};
for rr = 1:length(ROIs)
    ROI = ROIs{rr};
    figure;
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        % pRF
        prfpol = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.copol.prfs.nii.gz']));
        prfsig = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.cosig1.prfs.nii.gz']));
        prfco = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.co.prfs.nii.gz']));
        prfpeak = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.copeakt.prfs.nii.gz']));
        % CF
        cfpeak = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.cotrgpeakt.cfs.nii.gz']));
        cfshift = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.coshiftt.cfs.nii.gz']));
        if strcmp(srcROI,'volume')
            mask = load_nifti(fullfile(session_dir,[hemi '.' ROI '.mask.nii.gz']));
            maskind = (mask.vol > 0);
        elseif strcmp(srcROI,'cortex');
            mask = load_nifti(fullfile(session_dir,[hemi '.areas_pRF.nii.gz']));
            if strcmp(ROI,'V1')
                maskind = ((mask.vol == -1 | mask.vol == 1) & prfco.vol > 0.3873 & prfsig.vol > 0.25);
                dind = (mask.vol == -1 & prfco.vol > 0.3873 & prfsig.vol > 0.25);
                vind = (mask.vol == 1 & prfco.vol > 0.3873 & prfsig.vol > 0.25);
            elseif strcmp(ROI,'V2')
                maskind = ((mask.vol == -2 | mask.vol == 2) & prfco.vol > 0.3873 & prfsig.vol > 0.25);
                dind = (mask.vol == -2 & prfco.vol > 0.3873 & prfsig.vol > 0.25);
                vind = (mask.vol == 2 & prfco.vol > 0.3873 & prfsig.vol > 0.25);
                if strcmp(hemi,'lh')
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
                if strcmp(hemi,'lh')
                    dpolind = prfpol.vol > 0 & prfpol.vol < pi/2;
                    vpolind = prfpol.vol < 0 & prfpol.vol > -pi/2;
                elseif strcmp(hemi,'rh')
                    dpolind = prfpol.vol > pi/2 & prfpol.vol < pi;
                    vpolind = prfpol.vol < -pi/2 & prfpol.vol > -pi;
                end
                dind(~dpolind) = 0;
                vind(~vpolind) = 0;
            end
        end
        hold on;
        %
        if strcmp(srcROI,'volume')
            % ecc
            subplot(1,2,1 + (hh-1));
            difft = cfpeak.vol(maskind) - prfpeak.vol(maskind);
            plot(difft,cfshift.vol(maskind),'.');hold on;
            xlim([-5 5]);ylim([-5 5]);
            xlabel('HRF peak difference (seconds)','FontSize',10);
            ylabel('CF temporal shift (seconds)','FontSize',10);
            axis square
            title(hemi);
            annotation('textbox', tboxcoord{1 + 2*(hh-1)},...
                'String', ['mean HRF peak difference = ' num2str(mean(difft))]);
            annotation('textbox', tboxcoord{2 + 2*(hh-1)},...
                'String', ['mean CF temporal shift = ' num2str(mean(cfshift.vol(maskind)))]);
        elseif strcmp(srcROI,'cortex');
            % ecc dorsal
            subplot(2,2,1 + (hh-1));
            difft = cfpeak.vol(dind) - prfpeak.vol(dind);
            plot(difft,cfshift.vol(dind),[colors{hh + hh-1} '.']);hold on;
            xlim([-5 5]);ylim([-5 5]);
            xlabel('HRF peak difference (seconds)','FontSize',10);
            ylabel('CF temporal shift (seconds)','FontSize',10);
            axis square
            title([hemi ' dorsal']);
            annotation('textbox', tboxcoord{1 + 4*(hh-1)},...
                'String', ['mean HRF peak difference = ' num2str(mean(difft))]);
            annotation('textbox', tboxcoord{2 + 4*(hh-1)},...
                'String', ['mean CF temporal shift = ' num2str(mean(cfshift.vol(dind)))]);
            % ecc ventral
            subplot(2,2,3 + (hh-1));
            difft = cfpeak.vol(vind) - prfpeak.vol(vind);
            plot(difft,cfshift.vol(vind),[colors{hh + hh} '.']);hold on;
            xlim([-5 5]);ylim([-5 5]);
            xlabel('HRF peak difference (seconds)','FontSize',10);
            ylabel('CF temporal shift (seconds)','FontSize',10);
            axis square
            title([hemi ' ventral']);
            annotation('textbox', tboxcoord{3 + 4*(hh-1)},...
                'String', ['mean HRF peak difference = ' num2str(mean(difft))]);
            annotation('textbox', tboxcoord{4 + 4*(hh-1)},...
                'String', ['mean CF temporal shift = ' num2str(mean(cfshift.vol(vind)))]);
        end
    end
end
