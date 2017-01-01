function plot_CF_visual_field(session_dir,srcROI,trgROI,map_type,ROIs)

% Creates scatter plots of CF data. 
%
%   Usage:
%   plot_pRF_visual_field(session_dir,srcROI,ROIs)
%
%   Examples:
%   plot_pRF_visual_field(session_dir,'cortex',{'V1'});
%   plot_pRF_visual_field(session_dir,'volume');
%   
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
colors = {[1 0 0] [0 1 0] [0 0 1] [0 1 1]};
figure;
ct = 0;
for rr = 1:length(ROIs)
    ct = ct + 1;
    cct = 0;
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        eccfile = [hemi '.' srcROI '.' trgROI '.' map_type '.avg.coecc.cfs.nii.gz'];
        polfile = [hemi '.' srcROI '.' trgROI '.' map_type '.avg.copol.cfs.nii.gz'];
        sigfile = [hemi '.' srcROI '.' trgROI '.' map_type '.avg.cosig1.cfs.nii.gz'];
        cofile = [hemi '.' srcROI '.' trgROI '.' map_type '.avg.co.cfs.nii.gz'];
        ecc = load_nifti(fullfile(session_dir,eccfile));
        pol = load_nifti(fullfile(session_dir,polfile));
        sig = load_nifti(fullfile(session_dir,sigfile));
        co = load_nifti(fullfile(session_dir,cofile));
        [x,y] = pol2cart(pol.vol,ecc.vol);
        y = -y; % flip conversion
        sig = sig.vol;
        co = co.vol;
        cct = cct + 1;
        subplot(1,2,ct);
        ROI = ROIs{rr};
        if strcmp(srcROI,'volume')
            mask = load_nifti(fullfile(session_dir,[hemi '.' ROI '.mask.nii.gz']));
            maskind = (mask.vol > 0);
        elseif strcmp(srcROI,'cortex');
            mask = load_nifti(fullfile(session_dir,[hemi '.areas_pRF.nii.gz']));
            if strcmp(ROI,'V1')
                maskind = ((mask.vol == -1 | mask.vol == 1) & co > 0.3873 & sig > 0.25);
                dind = (mask.vol == -1 & co > 0.3873 & sig > 0.25);
                vind = (mask.vol == 1 & co > 0.3873 & sig > 0.25);
            elseif strcmp(ROI,'V2')
                maskind = ((mask.vol == -2 | mask.vol == 2) & co > 0.3873 & sig > 0.25);
                dind = (mask.vol == -2 & co > 0.3873 & sig > 0.25);
                vind = (mask.vol == 2 & co > 0.3873 & sig > 0.25);
            elseif strcmp(ROI,'V3')
                maskind = ((mask.vol == -3 | mask.vol == 3) & co > 0.3873 & sig > 0.25);
                dind = (mask.vol == -3 & co > 0.3873 & sig > 0.25);
                vind = (mask.vol == 3 & co > 0.3873 & sig > 0.25);
            end
        end
        % adjust sig
        xlim([-90 90]);
        ylim([-90 90]);
        axis square;
        currentunits = get(gca,'Units');
        set(gca, 'Units', 'Points');
        axpos = get(gca,'Position');
        set(gca, 'Units', currentunits);
        pts_per_angle = (axpos(3)^2)/diff(xlim);
        sigplot = sig*pts_per_angle; % Calculate sigma size in points^2
        hold on;
        %
        if strcmp(srcROI,'volume')
            scatter(x(maskind),y(maskind),sigplot(maskind),colors{cct});hold on;
        elseif strcmp(srcROI,'cortex');
            scatter(x(dind),y(dind),sigplot(dind),colors{cct + hh-1});hold on;
            scatter(x(vind),y(vind),sigplot(vind),colors{cct + hh});hold on;
        end
        grid on;
        title(ROI,'FontSize',30);
    end
    if strcmp(srcROI,'volume');
        legend(gca,'Left Hemisphere','Right Hemisphere');
    elseif strcmp(srcROI,'cortex');
        legend(gca,'Left Dorsal','Left Ventral', 'Right Dorsal', 'Right Ventral');
    end
    xlabel('Visual Angle (degrees)','FontSize',20);
    ylabel('Visual Angle (degrees)','FontSize',20);
end
