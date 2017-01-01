function [means,sems] = psc_cope(session_dir,subject_name,runNums,hemi,template,func,ROI,Fthresh,ctx_ROIs,copeNames,SUBJECTS_DIR)

% Gets the mean and SEM percent signal change for an ROI across runs, for 
%   each cope in a feat directory
%
%   Usage:
%   [means,sems] = psc_cope(session_dir,subject_name,hemi,template,func,ROI,Fthresh,ctx_ROIs,HzNames,SUBJECTS_DIR)
%
%   Written by Andrew S Bock Sep 2015

%% set defaults
if ~exist('template','var')
    template = 'anat';
end
if ~exist('func','var')
    func = 'brf.tf';
end
if ~exist('Fthresh','var')
    Fthresh = 4;
end
if ~exist('hemi','var')
    hemi = 'mh';
end
if ~exist('ctx_ROIs','var')
    ctx_ROIs = [5 15 50];
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
if ~exist('copeNames','var')
    copeNames = {...
        'Hz2' ...
        'Hz4' ...
        'Hz8' ...
        'Hz16' ...
        'Hz32' ...
        'Hz64' ...
        };
end
%% Find bold directories
b = find_bold(session_dir);
%% Get ROI indices
Fstatmask = load_nifti(fullfile(session_dir,[func '.Fstat.anat.nii.gz']));
switch template
    case 'anat'
        eccname = 'ecc';
        areaname = 'areas';
    case 'pRF'
        eccname = 'ecc_pRF';
        areaname = 'areas_pRF';
end
switch ROI
    case 'SC'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.SC.nii.gz']));
        ROIind = find(ROImask.vol>0 & Fstatmask.vol > Fthresh);
    case 'LGN'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.dots.LGN.nii.gz']));
        ROIind = find(ROImask.vol>0 & Fstatmask.vol > Fthresh);
    case 'MT'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.MT.nii.gz']));
        ROIind = find(ROImask.vol>0 & Fstatmask.vol > Fthresh);
    case 'pulvinar'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.pulvinar.nii.gz']));
        ROIind = find(ROImask.vol>0 & Fstatmask.vol > Fthresh);
    case 'V1low'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.' eccname '.vol.nii.gz']));
        areasmask = load_nifti(fullfile(session_dir,[hemi '.' areaname '.vol.nii.gz']));
        ROIind = find(ROImask.vol<=ctx_ROIs(1) & (areasmask.vol==1 | areasmask.vol==-1) ...
            & ~isnan(ROImask.vol) & Fstatmask.vol > Fthresh);
    case 'V1mid'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.' eccname '.vol.nii.gz']));
        areasmask = load_nifti(fullfile(session_dir,[hemi '.' areaname '.vol.nii.gz']));
        ROIind = find(ROImask.vol>ctx_ROIs(1) & ROImask.vol<=ctx_ROIs(2) & ...
            (areasmask.vol==1 | areasmask.vol==-1) & ~isnan(ROImask.vol) & Fstatmask.vol > Fthresh);
    case 'V1high'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.' eccname '.vol.nii.gz']));
        areasmask = load_nifti(fullfile(session_dir,[hemi '.' areaname '.vol.nii.gz']));
        ROIind = find(ROImask.vol>ctx_ROIs(2) & ROImask.vol<=ctx_ROIs(3) & ...
            (areasmask.vol==1 | areasmask.vol==-1) & ~isnan(ROImask.vol) & Fstatmask.vol > Fthresh);
    case 'V2low'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.' eccname '.vol.nii.gz']));
        areasmask = load_nifti(fullfile(session_dir,[hemi '.' areaname '.vol.nii.gz']));
        ROIind = find(ROImask.vol<=ctx_ROIs(1) & (areasmask.vol==2 | areasmask.vol==-2) & ...
            ~isnan(ROImask.vol) & Fstatmask.vol > Fthresh);
    case 'V2mid'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.' eccname '.vol.nii.gz']));
        areasmask = load_nifti(fullfile(session_dir,[hemi '.' areaname '.vol.nii.gz']));
        ROIind = find(ROImask.vol>ctx_ROIs(1) & ROImask.vol<=ctx_ROIs(2) & ...
            (areasmask.vol==2 | areasmask.vol==-2) & ~isnan(ROImask.vol) & Fstatmask.vol > Fthresh);
    case 'V2high'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.' eccname '.vol.nii.gz']));
        areasmask = load_nifti(fullfile(session_dir,[hemi '.' areaname '.vol.nii.gz']));
        ROIind = find(ROImask.vol>ctx_ROIs(2) & ROImask.vol<=ctx_ROIs(3) & ...
            (areasmask.vol==2 | areasmask.vol==-2) & ~isnan(ROImask.vol) & Fstatmask.vol > Fthresh);
    case 'V3low'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.' eccname '.vol.nii.gz']));
        areasmask = load_nifti(fullfile(session_dir,[hemi '.' areaname '.vol.nii.gz']));
        ROIind = find(ROImask.vol<=ctx_ROIs(1) & (areasmask.vol==3 | areasmask.vol==-3) & ...
            ~isnan(ROImask.vol) & Fstatmask.vol > Fthresh);
    case 'V3mid'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.' eccname '.vol.nii.gz']));
        areasmask = load_nifti(fullfile(session_dir,[hemi '.' areaname '.vol.nii.gz']));
        ROIind = find(ROImask.vol>ctx_ROIs(1) & ROImask.vol<=ctx_ROIs(2) & ...
            (areasmask.vol==3 | areasmask.vol==-3) & ~isnan(ROImask.vol) & Fstatmask.vol > Fthresh);
    case 'V3high'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.' eccname '.vol.nii.gz']));
        areasmask = load_nifti(fullfile(session_dir,[hemi '.' areaname '.vol.nii.gz']));
        ROIind = find(ROImask.vol>ctx_ROIs(2) & ROImask.vol<=ctx_ROIs(3) & ...
            (areasmask.vol==3 | areasmask.vol==-3) & ~isnan(ROImask.vol) & Fstatmask.vol > Fthresh);
    case 'V2/V3low'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.' eccname '.vol.nii.gz']));
        areasmask = load_nifti(fullfile(session_dir,[hemi '.' areaname '.vol.nii.gz']));
        ROIind = find(ROImask.vol<=ctx_ROIs(1) & (areasmask.vol==2 | areasmask.vol==-2 ...
            | areasmask.vol==3 | areasmask.vol==-3) & ~isnan(ROImask.vol) & Fstatmask.vol > Fthresh);
    case 'V2/V3mid'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.' eccname '.vol.nii.gz']));
        areasmask = load_nifti(fullfile(session_dir,[hemi '.' areaname '.vol.nii.gz']));
        ROIind = find(ROImask.vol>ctx_ROIs(1) & ROImask.vol<=ctx_ROIs(2) & ...
            (areasmask.vol==2 | areasmask.vol==-2 | areasmask.vol==3 | areasmask.vol==-3) ...
            & ~isnan(ROImask.vol) & Fstatmask.vol > Fthresh);
    case 'V2/V3high'
        ROImask = load_nifti(fullfile(session_dir,[hemi '.' eccname '.vol.nii.gz']));
        areasmask = load_nifti(fullfile(session_dir,[hemi '.' areaname '.vol.nii.gz']));
        ROIind = find(ROImask.vol>ctx_ROIs(2) & ROImask.vol<=ctx_ROIs(3) & ...
            (areasmask.vol==2 | areasmask.vol==-2 | areasmask.vol==3 | areasmask.vol==-3) ...
            & Fstatmask.vol > Fthresh);
end
%% Pull out beta values
ct = 0;
for i = runNums
    ct = ct + 1;
    for j = 1:length(copeNames)
        copefile = fullfile(session_dir,b{i},[func '.feat'],'stats',['cope' num2str(j) '.nii.gz']);
        copeout = fullfile(session_dir,b{i},[func '.feat'],'stats',['cope' num2str(j) '.anat.nii.gz']);
        if ~exist(copeout,'file');
            targ_vol = fullfile(SUBJECTS_DIR,subject_name,'mri','T1.mgz');
            tmpreg = listdir(fullfile(session_dir,b{i},'*bbreg.dat'),'files');
            bbreg_out_file = fullfile(session_dir,b{i},tmpreg{1}); % name registration file
            [~,~] = system(['mri_vol2vol --mov ' copefile ...
                ' --targ ' targ_vol ' --reg ' bbreg_out_file ...
                ' --o ' copeout ' --nearest']);
        end
        meanfile = fullfile(session_dir,b{i},[func '.feat'],'mean_func.nii.gz');
        meanout = fullfile(session_dir,b{i},[func '.feat'],'mean_func.anat.nii.gz');
        if ~exist(meanout,'file');
            targ_vol = fullfile(SUBJECTS_DIR,subject_name,'mri','T1.mgz');
            tmpreg = listdir(fullfile(session_dir,b{i},'*bbreg.dat'),'files');
            bbreg_out_file = fullfile(session_dir,b{i},tmpreg{1}); % name registration file
            [~,~] = system(['mri_vol2vol --mov ' meanfile ...
                ' --targ ' targ_vol ' --reg ' bbreg_out_file ...
                ' --o ' meanout ' --nearest']);
        end
        ctmp = load_nifti(copeout);
        mtmp = load_nifti(meanout);
        psctmp = (ctmp.vol./mtmp.vol)*100; % convert to percent signal change
        psctmp(psctmp==inf | psctmp==-inf) = nan; % set inf/-inf to nan
        eval([copeNames{j} '(ct) = nanmean(psctmp(ROIind));']);
    end
end
%% Create TTF
for j = 1:length(copeNames)
    eval(['means(j) = nanmean(' copeNames{j} ');']);
    eval(['sems(j) = nanstd(' copeNames{j} ') / sqrt(length(' copeNames{j} '));']);
end