function FIR_assemble(session_dir, subject_name, subj_name, output_dir, copeNames, runNums, hemis, ROIs, funcs, condition, startingCope)
% FIR_assemble(session_dir, subject_name, subj_name, output_dir, copeNames, ...
% runNums, hemis, ROIs, funcs, condition, startingCope)
%
% Saves out FIR responses and means
%
% <GF> Please add input arguments and usage here
%
% Input arguments:
% ================
%
%   session_dir : path to the session directory;
%   subject_name : Freesurfer subject name;
%   subj_name : subject name (will appear in the title of the figures);
%   output_dir : path to the output directory;
%   copeNames : names of the copes;
%   runNums : number of bold runs to include in the analysis (must be in the same
%   session_dir);
%   hemis : hemispheres to include in the analysis;
%   ROIs : Regions Of Interest to include in the analysis; 
%   funcs : smoothed functional data to include in the analysis (e.g.
%   wdrf.tf);
%   condition : condition to include in the analysis;
%   startingCope : from which cope the analysis starts (if not specified,
%   the analysis will start from the first cope).
%
% Usage:
% ======
%   session_dir = '/data/jag/MELA/MelanopsinMR/HERO_xxx1/060616'';
%   subject_name = 'HERO_xxx1_MaxMel';
%   subj_name = 'HERO_xxx1';
%   output_dir = '/data/jag/MELA/MelanopsinMR/Results';
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
%   runNums = 1:9;
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
%   condition = 'MelPulses_400pct';
%   startingCope = 15;
%   FIR_assemble(session_dir, subject_name, subj_name, output_dir, copeNames, ...
%   runNums, hemis, ROIs, funcs, condition, startingCope) 
%
%
% 7/2/2016  gf, ms      Written and commented.

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
            firFig = FIR_plot(mean,sem,ROI,condition,hemi,funcName, subj_name, runNums);
            if ~exist (fullfile(session_dir, 'FIR_figures'),'dir')
                mkdir (session_dir, 'FIR_figures');
            end
            if ~exist (fullfile(output_dir, 'FIR_figures'),'dir')
                mkdir (output_dir, 'FIR_figures');
            end
            savefig(fullfile(session_dir, 'FIR_figures', [subj_name '_' condition '_' hemi '_' ROI '_' func '.fig'])); %save .fig on cluster
            
            % Adjust the figure
            adjustPlot(firFig);
            saveas(firFig, fullfile(output_dir,'FIR_figures', [subj_name '_' condition '_' hemi '_' ROI '_' func '.pdf']), 'pdf');%save .pdf on dropbox
            close(firFig);
            
            % save means
            if ~exist (fullfile(session_dir, 'CSV_datafiles'),'dir')
                mkdir (session_dir, 'CSV_datafiles');
            end
            if ~exist (fullfile(output_dir, 'CSV_datafiles'),'dir')
                mkdir (output_dir, 'CSV_datafiles');
            end
            fileNameM = [subj_name '_' condition '_' hemi '_' ROI '_' func '_mean.csv'];
            fileNameS = [subj_name '_' condition '_' hemi '_' ROI '_' func '_SEM.csv'];
            csvwrite(fullfile(output_dir, 'CSV_datafiles', fileNameM), mean);
            csvwrite(fullfile(output_dir, 'CSV_datafiles', fileNameS), sem);
            csvwrite(fullfile(session_dir, 'CSV_datafiles', fileNameM), mean);
            csvwrite(fullfile(session_dir, 'CSV_datafiles', fileNameS), sem);
        end
    end
end