function project_copes(session_dir, subject_name, copeNames, hemis, ROIs, funcs)
% project_copes(session_dir, subject_name, copeNames, hemis, ROIs, funcs)
% 
%  Loop wrapper function for psc_copes.
% 
% Input arguments:
% ================
%
%   session_dir : path to the session directory; 
%   subject_name : Freesurfer subject name;
%   copeNames : names of the copes;
%   hemis : hemispheres to include in the analysis;
%   ROIs : Regions Of Interest to include in the analysis; 
%   funcs : smoothed functional data to include in the analysis (e.g.
%   wdrf.tf);
%
% Usage:
% ======
%   session_dir = '/data/jag/MELA/MelanopsinMR/HERO_xxx1/060616'';
%   subject_name = 'HERO_xxx1_MaxMel';
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
%   project_copes(session_dir, subject_name, copeNames, hemis, ROIs, funcs)
%
%
% 7/2/2016  gf, ms      Written and commented.

% Iterate over the hemis and ROIs
for hh = 1:length(hemis)
    for jj = 1:length(ROIs)
        ROI = ROIs{jj};
        % Get ROIind
        areas = load_nifti(fullfile(session_dir,[hemi '.areas.vol.nii.gz'])); % both hemis
        ecc = load_nifti(fullfile(session_dir,[hemi '.ecc.vol.nii.gz'])); % both hemis
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
        
        % Projects copes
        for ff = 1:length(funcs)
            func = funcs{ff};
            psc_cope(session_dir,subject_name,runNums,func,ROIind, copeNames);
        end
    end
end
