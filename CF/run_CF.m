function run_CF(session_dir,subject_name,runNum,hemi,srcROI,trgROI,func,DoG,seedSig1,seedSig2,seedSig3,seedSig4)

% Run the connetive field pipeline
%
%   Usage:
%   run_CF(session_dir,subject_name,runNum,hemi,srcROI,trgROI,func,DoG,seedSig1,seedSig2,seedSig3,seedSig4)
%
%   Written by Andrew S Bock May 2015

%% Find bold run directories
d = find_bold(session_dir);
%% Set up defaults
if ~exist('runNum','var')
    runNum = 1:length(d);
end
if ~exist('hemi','var')
    hemi = 'lh';
end
if ~exist('srcROI','var')
    srcROI = 'volume';
end
if ~exist('trgROI','var')
    trgROI='prf_V1';
end
if ~exist('DoG','var');
    DoG = 1;
end
if ~exist('seedSig1','var');
    %seedSig1 = (.5:.5:10)';
    %seedSig1 = (0.5:.5:5)';
    %seedSig1 = (1:1:10)';
    seedSig1 = (0.5:0.5:15)';
end
if ~exist('seedSig2','var');
    %seedSig2 = (1:.5:3)';
    %seedSig2 = (1:.5:3)';
    %seedSig2 = (1.25:.5:3)';
    seedSig2 = 0;
end
if ~exist('seedSig3','var');
    %seedSig3 = (0:0.125:1)';
    %seedSig3 = (0:0.25:1)';
    seedSig3 = 1;
end
if ~exist('seedSig4','var');
    %seedSig4 = (0:0.125:1)';
    %seedSig4 = (0:0.25:1)';
    seedSig4 = 0;
end
seedSig = {seedSig1 seedSig2 seedSig3 seedSig4};

%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,subject_name,runNum,hemi,srcROI,trgROI,DoG,seedSig1,seedSig2,seedSig3,seedSig4);

%% Run the CF pipeline
for rr = runNum;
    disp(['Generating ' srcROI ' CFs for ' d{rr}]);
    cd(fullfile(session_dir,d{rr}));
    % Get source indices
    switch trgROI
        case 'V1'
            areas = load_nifti(fullfile(session_dir,[hemi '.areas.nii.gz']));
            ecc = load_nifti(fullfile(session_dir,[hemi '.ecc.nii.gz']));
        case 'prf_V1'
            areas = load_nifti(fullfile(session_dir,[hemi '.areas_pRF.nii.gz']));
            ecc = load_nifti(fullfile(session_dir,[hemi '.ecc_pRF.nii.gz']));
    end
    switch srcROI
        case 'V3'
            V1ind = find(areas.vol<=1 & areas.vol >=-1);
            V1_3ind = find(areas.vol<=3 & areas.vol >=-3);
            V1_3ind(V1ind) = 0;
            srcind = find(V1_3ind);
        case 'cortex'
            srcind = 1:length(areas.vol); % entire cortex
        case 'volume'
            %binfile = fullfile(session_dir,d{rr},'brf.aseg.gm.nii.gz');
            binfile = fullfile(session_dir,d{rr},'single_TR.nii.gz');
            src = load_nifti(binfile);
            srcind = find(src.vol > 0 & ~isnan(src.vol));
        case 'LGN'
            binfile = fullfile(session_dir,d{rr},[hemi '.LGN.nii.gz']);
            src = load_nifti(binfile);
            srcind = find(src.vol > 0 & ~isnan(src.vol));
    end
    % Get target indices
    trgind = find(areas.vol<=1 & areas.vol >=-1);
    if strcmp(srcROI,'cortex')
        srcfile = fullfile(session_dir,d{rr},[func '_surf.' hemi '.nii.gz']);
    else
        srcfile = fullfile(session_dir,d{rr},[func '.nii.gz']);
    end
    trgfile = fullfile(session_dir,d{rr},['s5.dbrf.tf_surf.' hemi '.nii.gz']);
    do_CF(session_dir,subject_name,rr,srcfile,trgfile,srcind,trgind,srcROI,trgROI,seedSig,hemi,DoG);
end
