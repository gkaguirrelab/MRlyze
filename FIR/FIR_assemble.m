function FIR_assemble(session_dir, subject_name, subj_name, output_dir, SUBJECTS_DIR, copeNames, hemis, ROIs, funcs, condition)
% Saves out FIR responses and means


% FIR_assemble(session_dir, subject_name, subj_name, output_dir, SUBJECTS_DIR, copeNames, hemis, ROIs, funcs, funcNames, Conditions, Runs)


for hh = 1:length(hemis)
    hemi = hemis{hh};
    for jj = 1:length(ROIs)
        ROI = ROIs{jj};
        % Get ROIind
        areas = load_nifti(fullfile(session_dir,[hemi '.areas.vol.nii.gz'])); 
        ecc = load_nifti(fullfile(session_dir,[hemi '.ecc.vol.nii.gz'])); 
        switch ROI
            case 'V1'
                ROIind = find(abs(areas.vol)==1 & (ecc.vol>5 & ecc.vol<=30));
            case 'V2andV3'
                ROIind = find(abs(areas.vol)==2 | abs(areas.vol)==3);
            case 'LGN'
                if strcmp(hemi,'mh')
                    lgn_rh = load_nifti(fullfile(session_dir, 'rh.LGN.prob.nii.gz'));
                    lgn_lh = load_nifti(fullfile(session_dir, 'lh.LGN.prob.nii.gz'));
                    ROIind = find((lgn_rh.vol)>=25 | (lgn_lh.vol)>=25);
                else
                    lgn = load_nifti(fullfile(session_dir, [hemi '.LGN.prob.nii.gz']));
                    ROIind = find((lgn.vol)>=25);
                end
        end
        
        %% Get means
        for ff = 1:length(funcs)
            func = funcs{ff};
            funcName = funcs{ff};
            condName = condition;
            [means{i},sems{i}] = psc_cope_get_means(session_dir,subject_name,runNums,func,ROIind, copeNames);
            % Plot and save figures
            FIR_plot(means{i},sems{i},ROI,condName,hemi,funcName, subj_name, runNums);
            if ~exist (fullfile(session_dir, 'FIR_figures'),'dir')
                mkdir (session_dir, 'FIR_figures');
            end
            if ~exist (fullfile(output_dir, 'FIR_figures'),'dir')
                mkdir (output_dir, 'FIR_figures');
            end
            savefig(fullfile(session_dir, 'FIR_figures', [subj_name '_' condName '_' hemi '_' ROI '_' func '.fig'])); %save .fig on cluster
            set(gcf, 'PaperPosition', [0 0 7 7]);
            set(gcf, 'PaperSize', [7 7]);
            saveas(gcf, fullfile(output_dir,'FIR_figures', [subj_name '_' condName '_' hemi '_' ROI '_' func]), 'pdf');%save .pdf on dropbox
            close all;
            % save means
            if ~exist (fullfile(session_dir, 'CSV_datafiles'),'dir')
                mkdir (session_dir, 'CSV_datafiles');
            end
            if ~exist (fullfile(output_dir, 'CSV_datafiles'),'dir')
                mkdir (output_dir, 'CSV_datafiles');
            end
            fileNameM = [subj_name '_' condName '_' hemi '_' ROI '_' func '_' 'mean.csv'];
            fileNameS = [subj_name '_' condName '_' hemi '_' ROI '_' func '_' 'SEM.csv'];
            csvwrite ((fullfile(output_dir,'CSV_datafiles', fileNameM)), means);
            csvwrite ((fullfile(output_dir,'CSV_datafiles', fileNameS)), sems);
            csvwrite ((fullfile(session_dir,'CSV_datafiles', fileNameM)), means);
            csvwrite ((fullfile(session_dir,'CSV_datafiles', fileNameS)), sems);
            
        end
    end
end