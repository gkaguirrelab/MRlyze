function project_copes(session_dir, subject_name, SUBJECTS_DIR, copeNames, hemis, ROIs, funcs, funcNames, Conditions, Runs, PROJECT_TEMPLATE)

copeNames = {...
    'Sec00' ...
    'Sec01' ...
    'Sec02' ...
    'Sec03' ...
    'Sec04' ...
    'Sec05' ...
    'Sec06' ...
    'Sec07' ...
    'Sec08' ...
    'Sec09' ...
    'Sec10' ...
    'Sec11' ...
    'Sec12' ...
    'Sec13' ...
    };

%%
for hh = 1:length(hemis)
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
        
        % Projects copes
        for ff = 1:length(funcs)
            func = funcs{ff};
            for i = 1:length(Conditions)
                runNums = Runs{i};
                psc_cope(session_dir,subject_name,runNums,func,ROIind, copeNames);
            end
        end
    end
end