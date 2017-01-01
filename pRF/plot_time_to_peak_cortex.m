function plot_time_to_peak_cortex(session_dir)

srcROI = 'cortex';
trgROI = 'prf_V1';
map_type = 'movie';
hemis = {'lh' 'rh'};
ROIs = {'V2' 'V3'};

for rr = 1:length(ROIs)
    figure(1);
    for hh = 1:length(hemis)
        % define variable names
        hemi = hemis{hh};
        ROI = ROIs{rr};
        mask_file = fullfile(session_dir,[hemi '.areas_pRF.nii.gz']);
        pRF_srcpeak = fullfile(session_dir,[hemi '.cortex.avg.copeakt.prfs.nii.gz']);
        pRF_trgpeak = fullfile(session_dir,[hemi '.cortex.avg.copeakt.prfs.nii.gz']);
        pRF_ecc = fullfile(session_dir,[hemi '.ecc_pRF.nii.gz']);
        pRF_pol = fullfile(session_dir,[hemi '.pol_pRF.nii.gz']);
        pRF_areas = fullfile(session_dir,[hemi '.areas_pRF.nii.gz']);
        CF_srcpeak = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.copeakt.cfs.nii.gz']);
        CF_ecc = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.coecc.cfs.nii.gz']);
        CF_pol = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.copol.cfs.nii.gz']);
        % load files
        mask = load_nifti(mask_file);
        pRFsrcpeak = load_nifti(pRF_srcpeak);
        pRFtrgpeak = load_nifti(pRF_trgpeak);
        pRFecc = load_nifti(pRF_ecc);
        pRFpol = load_nifti(pRF_pol);
        pRFareas = load_nifti(pRF_areas);
        CFsrcpeak = load_nifti(CF_srcpeak);
        CFecc = load_nifti(CF_ecc);
        CFpol = load_nifti(CF_pol);
        % Find CF target time to peak
        if strcmp(ROI,'V2')
            maskind = find(mask.vol==2 | mask.vol == -2);
        elseif strcmp(ROI,'V3')
            maskind = find(mask.vol==3 | mask.vol ==-3);
        end
        V1ind = find(pRFareas.vol >=-1 & pRFareas.vol<=1);
        V1ecc = pRFecc.vol(V1ind);
        V1pol = pRFpol.vol(V1ind);
        CFtrgpeak = CFsrcpeak;
        CFtrgpeak.vol = nan(size(CFtrgpeak.vol));
        pRFdiffpeak = CFtrgpeak;
        for mm = 1:length(maskind);
            tmpecc = CFecc.vol(maskind(mm));
            tmppol = CFpol.vol(maskind(mm));
            [~,eccdiffind] = sort(abs(tmpecc - V1ecc));
            [~,poldiffind] = sort(abs(tmppol - V1pol));
            [~,eccorder] = sort(eccdiffind);
            [~,polorder] = sort(poldiffind);
            sumorder = eccorder + polorder;
            [~,bestind] = min(sumorder);
            CFtrgpeak.vol(maskind(mm)) = pRFtrgpeak.vol(bestind);
        end
        
        pRFdiffpeak.vol(maskind) =  CFtrgpeak.vol(maskind) - pRFsrcpeak.vol(maskind);
         % new figure for difference values
        if rr ==1
            figure(1)
        else
            figure(3)
        end
        % Difference between V1 time to peak and ROI time to peak
        X = -4:0.5:4; % set X axis
        if hh == 1
            ct = 1;
        else
            ct = 2;
        end
        dN = histc(pRFdiffpeak.vol(maskind),X);
        subplot(2,2,ct);
        bar(X,dN,1)
        title(['Difference between V1 and ' ROI ' HRF time to peak: ' hemi],'FontSize',25);
        ylim([0 2000]);
        set(gca,'xlim',[min(X) max(X)],'XTick',X);
        axis square
        xlabel('Time (seconds)','FontSize',25)
        ylabel('Number of voxels','FontSize',25)
        hold on
        dmean = mean(pRFdiffpeak.vol(maskind));
        maxY = ylim;maxY = maxY(2);
        dH = bar(dmean,maxY,0.05,'FaceColor',[1 0 0]);
        legend(dH,['Mean = ' num2str(dmean)]);
        % Shift in ROI time to peak
        if hh == 1
            ct = 3;
        else
            ct = 4;
        end
        rN = histc(CFsrcpeak.vol(maskind),X);
        subplot(2,2,ct);
        bar(X,rN,1)
        title(['Shift in CF timecourse: ' ROI ' ' hemi],'FontSize',25);
        ylim([0 4000]);
        set(gca,'xlim',[min(X) max(X)],'XTick',X);
        axis square
        xlabel('Time (seconds)','FontSize',25)
        ylabel('Number of voxels','FontSize',25)
        hold on
        rmean = mean(CFsrcpeak.vol(maskind));
        maxY = ylim;maxY = maxY(2);
        rH = bar(rmean,maxY,0.05,'FaceColor',[1 0 0]);
        legend(rH,['Mean = ' num2str(rmean)]);
        
        % new figure for raw values
        if rr ==1
            figure(2)
        else
            figure(4)
        end
        % raw time to peak values V1
        X = 0:0.5:10; % set X axis
        if hh == 1
            ct = 1;
        else
            ct = 2;
        end
        vN = histc(pRFtrgpeak.vol(V1ind),X);
        subplot(2,2,ct);
        bar(X,vN,1)
        title(['HRF time to peak: V1 ' hemi],'FontSize',25);
        ylim([0 2000]);
        set(gca,'xlim',[min(X) max(X)],'XTick',X);
        axis square
        xlabel('Time (seconds)','FontSize',25)
        ylabel('Number of voxels','FontSize',25)
        hold on
        dmean = mean(pRFtrgpeak.vol(V1ind));
        maxY = ylim;maxY = maxY(2);
        dH = bar(dmean,maxY,0.05,'FaceColor',[1 0 0]);
        legend(dH,['Mean = ' num2str(dmean)]);
        % raw time to peak values (ROI)
        if hh == 1
            ct = 3;
        else
            ct = 4;
        end
        sN = histc(pRFsrcpeak.vol(maskind),X);
        subplot(2,2,ct);
        bar(X,sN,1)
        title(['HRF time to peak: ' ROI ' ' hemi],'FontSize',25);
        ylim([0 2000]);
        set(gca,'xlim',[min(X) max(X)],'XTick',X);
        axis square
        xlabel('Time (seconds)','FontSize',25)
        ylabel('Number of voxels','FontSize',25)
        hold on
        dmean = mean(pRFsrcpeak.vol(maskind));
        maxY = ylim;maxY = maxY(2);
        dH = bar(dmean,maxY,0.05,'FaceColor',[1 0 0]);
        legend(dH,['Mean = ' num2str(dmean)]);
    end
end