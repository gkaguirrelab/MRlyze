function compare_CF_pRF(session_dir,srcROI,trgROI,map_type,ROIs)

% Creates scatter plots of CF data.
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
for rr = 1:length(ROIs)
    ROI = ROIs{rr};
    figure;
    ct = 0;
    for hh = 1:length(hemis)
        ct = ct + 1;
        hemi = hemis{hh};
        % pRF
        prfecc = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.coecc.prfs.nii.gz']));
        prfpol = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.copol.prfs.nii.gz']));
        prfsig = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.cosig1.prfs.nii.gz']));
        prfco = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.co.prfs.nii.gz']));
        % CF
        cfecc = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.coecc.cfs.nii.gz']));
        cfpol = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.copol.cfs.nii.gz']));
        cfsig = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.cosig1.cfs.nii.gz']));
        cfco = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.co.cfs.nii.gz']));
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
            subplot(4,2,1 + (ct-1));
            plot(prfecc.vol(maskind),cfecc.vol(maskind),'.');hold on;
            xlim([0 35]);ylim([0 35]);
            xlabel('pRF Eccentricty (degrees)','FontSize',10);
            ylabel('CF Eccentricty (degrees)','FontSize',10);
            axis square
            corrval = corr2(prfecc.vol(maskind),cfecc.vol(maskind));
            title(hemi);
            legend(['corr = ' num2str(corrval)]);
            % pol
            subplot(4,2,3 + (ct-1));
            plot(prfpol.vol(maskind),cfpol.vol(maskind),'.');hold on;
            xlim([-pi pi]);ylim([-pi pi]);
            xlabel('pRF Polar Angle (radians)','FontSize',10);
            ylabel('CF Polar Angle (radians)','FontSize',10);
            axis square
            [corrval] = circ_corrcc(prfpol.vol(maskind),cfpol.vol(maskind));
            title(hemi);
            legend(['corr = ' num2str(corrval)]);
            % sig
            subplot(4,2,5 + (ct-1));
            plot(prfsig.vol(maskind),cfsig.vol(maskind),'.');hold on;
            xlim([0 10]);ylim([0 10]);
            xlabel('pRF Sigma (degrees)','FontSize',10);
            ylabel('CF Sigma (mm)','FontSize',10);
            axis square
            corrval = corr2(prfecc.vol(maskind),cfecc.vol(maskind));
            title(hemi);
            legend(['corr = ' num2str(corrval)]);
            % co
            subplot(4,2,7 + (ct-1));
            plot(prfco.vol(maskind),cfco.vol(maskind),'.');hold on;
            xlim([0 1]);ylim([0 1]);
            xlabel('pRF Correlation (fisher-z))','FontSize',10);
            ylabel('CF Correlation (fisher-z)','FontSize',10);
            axis square
            corrval = corr2(prfecc.vol(maskind),cfecc.vol(maskind));
            title(hemi);
            legend(['corr = ' num2str(corrval)]);
        elseif strcmp(srcROI,'cortex');
            % ecc dorsal
            subplot(4,4,1 + 2*(ct-1));
            plot(prfecc.vol(dind),cfecc.vol(dind),[colors{ct + hh-1} '.']);hold on;
            xlim([0 35]);ylim([0 35]);
            xlabel('pRF Eccentricty (degrees)','FontSize',10);
            ylabel('CF Eccentricty (degrees)','FontSize',10);
            axis square
            corrval = corr2(prfecc.vol(dind),cfecc.vol(dind));
            title([hemi ' dorsal']);
            legend(['corr = ' num2str(corrval)]);
            % ecc ventral
            subplot(4,4,2 + 2*(ct-1));
            plot(prfecc.vol(vind),cfecc.vol(vind),[colors{ct + hh} '.']);hold on;
            xlim([0 35]);ylim([0 35]);
            xlabel('pRF Eccentricty (degrees)','FontSize',10);
            ylabel('CF Eccentricty (degrees)','FontSize',10);
            axis square
            corrval = corr2(prfecc.vol(vind),cfecc.vol(vind));
            title([hemi ' ventral']);
            legend(['corr = ' num2str(corrval)]);
            % pol dorsal
            subplot(4,4,5 + 2*(ct-1));
            plot(prfpol.vol(dind),cfpol.vol(dind),[colors{ct + hh-1} '.']);hold on;
            xlim([-pi pi]);ylim([-pi pi]);
            xlabel('pRF Polar Angle (radians)','FontSize',10);
            ylabel('CF Polar Angle (radians)','FontSize',10);
            axis square
            [corrval] = circ_corrcc(prfpol.vol(dind),cfpol.vol(dind));
            title([hemi ' dorsal']);
            legend(['corr = ' num2str(corrval)]);
            % pol ventral
            subplot(4,4,6 + 2*(ct-1));
            plot(prfpol.vol(vind),cfpol.vol(vind),[colors{ct + hh} '.']);hold on;
            xlim([-pi pi]);ylim([-pi pi]);
            xlabel('pRF Polar Angle (radians)','FontSize',10);
            ylabel('CF Polar Angle (radians)','FontSize',10);
            axis square
            [corrval] = circ_corrcc(prfpol.vol(vind),cfpol.vol(vind));
            title([hemi ' ventral']);
            legend(['corr = ' num2str(corrval)]);
            % sig dorsal
            subplot(4,4,9 + 2*(ct-1));
            plot(prfsig.vol(dind),cfsig.vol(dind),[colors{ct + hh-1} '.']);hold on;
            xlim([0 10]);ylim([0 10]);
            xlabel('pRF Sigma (degrees)','FontSize',10);
            ylabel('CF Sigma (mm)','FontSize',10);
            axis square
            corrval = corr2(prfsig.vol(dind),cfsig.vol(dind));
            title([hemi ' dorsal']);
            legend(['corr = ' num2str(corrval)]);
            % sig ventral
            subplot(4,4,10 + 2*(ct-1));
            plot(prfsig.vol(vind),cfsig.vol(vind),[colors{ct + hh} '.']);hold on;
            xlim([0 10]);ylim([0 10]);
            xlabel('pRF Sigma (degrees)','FontSize',10);
            ylabel('CF Sigma (mm)','FontSize',10);
            axis square
            corrval = corr2(prfsig.vol(vind),cfsig.vol(vind));
            title([hemi ' ventral']);
            legend(['corr = ' num2str(corrval)]);
            % co dorsal
            subplot(4,4,13 + 2*(ct-1));
            plot(prfco.vol(dind),cfco.vol(dind),[colors{ct + hh-1} '.']);hold on;
            xlim([0 1]);ylim([0 1]);
            xlabel('pRF Correlation (fisher-z))','FontSize',10);
            ylabel('CF Correlation (fisher-z)','FontSize',10);
            axis square
            corrval = corr2(prfco.vol(dind),cfco.vol(dind));
            title([hemi ' dorsal']);
            legend(['corr = ' num2str(corrval)]);
            % co ventral
            subplot(4,4,14 + 2*(ct-1));
            plot(prfco.vol(vind),cfco.vol(vind),[colors{ct + hh} '.']);hold on;
            xlim([0 1]);ylim([0 1]);
            xlabel('pRF Correlation (fisher-z))','FontSize',10);
            ylabel('CF Correlation (fisher-z)','FontSize',10);
            axis square
            corrval = corr2(prfco.vol(vind),cfco.vol(vind));
            title([hemi ' ventral']);
            legend(['corr = ' num2str(corrval)]);
        end
    end
end
