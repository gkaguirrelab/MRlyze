function plot_time_to_peak_volume(session_dir,map_type)

srcROI = 'volume';
trgROI = 'prf_V1';
if ~exist('map_type','var')
    map_type = 'movie';
end
hemis = {'lh' 'rh'};
ROIs = {'LGN' 'SC'};
figct = 0;
for rr = 1:length(ROIs)
    for hh = 1:length(hemis)
        % load files
        hemi = hemis{hh};
        ROI = ROIs{rr};
        mask_file = load_nifti(fullfile(session_dir,[hemi '.' ROI '.mask.nii.gz']));
        shift_file = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.coshiftt.cfs.nii.gz']));
        srcHRF_file = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.avg.copeakt.prfs.nii.gz']));
        trgHRF_file = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.cotrgpeakt.cfs.nii.gz']));
       
        

        
        shift_file = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.coshiftt.cfs.nii.gz']));
        diff_file = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.codiffpeakt.cfs.nii.gz']));
        diffshift_file = load_nifti(fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.' map_type '.avg.codiffshift.cfs.nii.gz']));
        % Find CF target time to peak
        maskind = find(mask_file.vol>0);    
        % new figure for comparison plots
        if rr ==1
            figure(5)
        else
            figure(6)
        end
        if hh == 1
            ct = 1;
        else
            ct = 2;
        end
        subplot(1,2,ct);
        plot(diff_file.vol(maskind),shift_file.vol(maskind),'.');hold on
        xlim([-5 5]);
        ylim([-5 5]);
        axis square;lsline
        title([ROI ' ' hemi],'FontSize',25);
        xlabel('HRF time to peak difference','FontSize',25)
        ylabel('Temporal offset parameter in CF data','FontSize',25)
    end
end

%%
difffoo = trgHRF_file.vol(maskind) - srcHRF_file.vol(maskind);
figure;plot(difffoo,shift_file.vol(maskind),'.');
xlim([-3 3]);
ylim([-3 3]);
lsline;

figure;hist(shift_file.vol(maskind),100);
figure;hist(difffoo,100);


%%
close all
figure;hist(srcHRF_file.vol(maskind),100)
figure;hist(trgHRF_file.vol(maskind),100)
figure;hist(diff_file.vol(maskind),100)
