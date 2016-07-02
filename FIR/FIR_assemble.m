function FIR_assemble(session_dir, subject_name, subj_name, output_dir, copeNames, runNums, hemis, ROIs, funcs, condition, startingCope)
% FIR_assemble(session_dir, subject_name, subj_name, output_dir, copeNames, runNums, hemis, ROIs, funcs, condition, startingCope)
%
% Saves out FIR responses and means

if ~exist('startingCope','var')
    startingCope = 1;
end

% Iterate over hemis and ROIs
for hh = 1:length(hemis)
    hemi = hemis{hh};
    for jj = 1:length(ROIs)
        ROI = ROIs{jj};
        % Get ROIind
        areas = load_nifti(fullfile(session_dir,'anat_templates',[hemi '.areas.anat.vol.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,'anat_templates',[hemi '.ecc.anat.vol.nii.gz']));
        lgn = load_nifti(fullfile(session_dir,'anat_templates',[hemi '.LGN.nii.gz']));
        switch ROI
            case 'V1'
                ROIind = find(abs(areas.vol)==1 & (ecc.vol>5 & ecc.vol<=30)); 
            case 'V2andV3'
                ROIind = find(abs(areas.vol)==2 | abs(areas.vol)==3);
            case 'LGN'
                ROIind = find((lgn.vol)>0);
        end
        %% Get means
        for ff = 1:length(funcs)
            func = funcs{ff};
            funcName = funcs{ff};
            [mean,sem] = psc_cope_get_means(session_dir,subject_name,runNums,func,ROIind,copeNames,startingCope);
            
            % Plot and save figures
            FIR_plot(mean,sem,ROI,condition,hemi,funcName, subj_name, runNums);
            if ~exist (fullfile(session_dir, 'FIR_figures'),'dir')
                mkdir (session_dir, 'FIR_figures');
            end
            if ~exist (fullfile(output_dir, 'FIR_figures'),'dir')
                mkdir (output_dir, 'FIR_figures');
            end
            savefig(fullfile(session_dir, 'FIR_figures', [subj_name '_' condition '_' hemi '_' ROI '_' func '.fig'])); %save .fig on cluster
            
            % Adjust the figure
            adjustPlot(gcf);
            saveas(gcf, fullfile(output_dir,'FIR_figures', [subj_name '_' condition '_' hemi '_' ROI '_' func '.pdf']), 'pdf');%save .pdf on dropbox
            close all;
            
            % save means
            if ~exist (fullfile(session_dir, 'CSV_datafiles'),'dir')
                mkdir (session_dir, 'CSV_datafiles');
            end
            if ~exist (fullfile(output_dir, 'CSV_datafiles'),'dir')
                mkdir (output_dir, 'CSV_datafiles');
            end
            fileNameM = [subj_name '_' condition '_' hemi '_' ROI '_' func '_mean.csv'];
            fileNameS = [subj_name '_' condition '_' hemi '_' ROI '_' func '_SEM.csv'];
            csvwrite ((fullfile(output_dir,'CSV_datafiles', fileNameM)), mean);
            csvwrite ((fullfile(output_dir,'CSV_datafiles', fileNameS)), sem);
            csvwrite ((fullfile(session_dir,'CSV_datafiles', fileNameM)), mean);
            csvwrite ((fullfile(session_dir,'CSV_datafiles', fileNameS)), sem);
            
        end
    end
end