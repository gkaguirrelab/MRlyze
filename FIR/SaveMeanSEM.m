function SaveMeanSEM(session_dir, subject_name, subj_name, dropbox_dir, SUBJECTS_DIR, copeNames, hemis, ROIs, funcs, funcNames, Conditions, Runs)
% SaveMeanSEM(session_dir, subject_name, subj_name, dropbox_dir, ...
% SUBJECTS_DIR, copeNames, hemis, ROIs, funcs, funcNames, Conditions, Runs)
% 
%
% Input arguments:
% ================
%
%   session_dir : path to the session directory;
%   subject_name : Freesurfer subject name;
%   subj_name : subject name (will appear in the filenames);
%   dropbox_dir : path to the dropbox output directory;
%   SUBJECTS_DIR : freesurfer subjects directory;
%   copeNames : names of the copes;
%   hemis : hemispheres to include in the analysis;
%   ROIs : Regions Of Interest to include in the analysis; 
%   funcs : smoothed functional data to include in the analysis (e.g.
%   wdrf.tf);
%   funcsNames : name of the funcs (in case it is different from the
%   funcs);
%   Conditions : conditions to include in the analysis;
%   Runs: runs to include in the analysis.
%
% Usage:
% ======
%   session_dir = '/data/jag/MELA/MelanopsinMR/HERO_xxx1/060616'';
%   subject_name = 'HERO_xxx1_MaxMel';
%   subj_name = 'HERO_xxx1';
%   drobpox_dir = '/Users/giulia/Aguirre-Brainard-Dropbox/MELA_analysis/';
%   SUBJECTS_DIR = getenv(SUBJECTS_DIR);
%   copeNames = {...
%             'Sec00' ...
%             'Sec01' ...
%             'Sec02' ...
%             'Sec03' ...
%             'Sec04' ...
%             'Sec05' ...
%             'Sec06' ...
%             'Sec07' ...
%             'Sec08' ...
%             'Sec09' ...
%             'Sec10' ...
%             'Sec11' ...
%             'Sec12' ...
%             'Sec13' ...
%             };
%    hemis = {...
%             'mh'...
%             'lh'...
%             'rh'...
%             };
%   ROIs = {...
%             'V1' ...
%             'V2andV3'...
%             'LGN'...
%             }; 
%   funcs = {...
%     'wdrf.tf' ... %raw data 
%     's5.wdrf.tf' ... %5mm smoothed data
%     };
%   Conditions = { ...
%           'MelPulses_400pct';
%           };
%   Runs = 1:9;
% SaveMeanSEM(session_dir, subject_name, subj_name, dropbox_dir, ...
% SUBJECTS_DIR, copeNames, hemis, ROIs, funcs, funcNames, Conditions, Runs)
% 
%   
%
% 7/2/2016  gf, ms      Written and commented.

% Iterate over hemis
for hh = 1:length(hemis)
    hemi = hemis{hh};
    %LGN variables
    in_vol = fullfile('~/data/' , [hemi '.LGN.prob.nii.gz']);
    out_vol =  fullfile(session_dir, [hemi '.LGN.prob.nii.gz']);
    ref_vol = fullfile (SUBJECTS_DIR , subject_name, '/mri/T1.mgz');
    
    % Iterate over ROIs
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
            csvwrite(fullfile(dropbox_dir,'CSV_datafiles', fileNameM), means);
            csvwrite(fullfile(dropbox_dir, 'CSV_datafiles', fileNameS), sems);
            csvwrite(fullfile(session_dir,'CSV_datafiles', fileNameM), means);
            csvwrite(fullfile(session_dir, 'CSV_datafiles', fileNameS), sems);
            
        end
        
    end
end