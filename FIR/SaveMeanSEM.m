function SaveMeanSEM(session_dir, subject_name, subj_name, dropbox_dir, SUBJECTS_DIR, copeNames, hemis, ROIs, funcs, funcNames, Conditions, Runs)
% FIR_assemble(session_dir, subject_name, subj_name, dropbox_dir, SUBJECTS_DIR, copeNames, hemis, ROIs, funcs, funcNames, Conditions, Runs)

for hh = 1:length(hemis)
    hemi = hemis{hh};
    %LGN variables
    in_vol = fullfile('~/data/' , [hemi '.LGN.prob.nii.gz']);
    out_vol =  fullfile(session_dir, [hemi '.LGN.prob.nii.gz']);
    ref_vol = fullfile (SUBJECTS_DIR , subject_name, '/mri/T1.mgz');

    %%
    
    for jj = 1:length(ROIs)
        ROI = ROIs{jj};
        % Get ROIind
        areas = load_nifti(fullfile(session_dir,[hemi '.areas.vol.nii.gz'])); % both hemis
        ecc = load_nifti(fullfile(session_dir,[hemi '.ecc.vol.nii.gz'])); % both hemis
        switch ROI
            case 'V1all'
                ROIind = find(abs(areas.vol)==1); % all of V1
            case 'V1low'
                ROIind = find(abs(areas.vol)==1 & ecc.vol<=5);
            case 'V1mid'
                ROIind = find(abs(areas.vol)==1 & (ecc.vol>5 & ecc.vol<=15));
            case 'V1high'
                ROIind = find(abs(areas.vol)==1 & (ecc.vol>15 & ecc.vol<=40));
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
            case 'V1'
                ROIind = find(abs(areas.vol)==1 & (ecc.vol>5 & ecc.vol<=30));
        end        
        
        %% Get means
        for ff = 1:length(funcs)
            func = funcs{ff};
            funcName = funcNames{ff};
            for i = 1:length(Conditions)
                condName = Conditions{i};
                runNums = Runs{i};
                fileNameM = [subj_name '_' condName '_' hemi '_' ROI '_' func '_' 'mean.csv'];
                fileNameS = [subj_name '_' condName '_' hemi '_' ROI '_' func '_' 'SEM.csv'];
                 if ~exist (fullfile(session_dir, 'CSV_datafiles'),'dir')
                    mkdir (session_dir, 'CSV_datafiles');
                end
                if ~exist (fullfile(dropbox_dir, 'CSV_datafiles'),'dir')
                    mkdir (dropbox_dir, 'CSV_datafiles');
                end
                
                [means{i},sems{i}] = psc_cope_get_means(session_dir,subject_name,runNums,func,ROIind, copeNames);  
            end
            csvwrite ((fullfile(dropbox_dir,'CSV_datafiles', fileNameM)), means);
            csvwrite ((fullfile(dropbox_dir, 'CSV_datafiles', fileNameS)), sems);
            csvwrite ((fullfile(session_dir,'CSV_datafiles', fileNameM)), means);
            csvwrite ((fullfile(session_dir, 'CSV_datafiles', fileNameS)), sems);
            
        end
        
    end
end