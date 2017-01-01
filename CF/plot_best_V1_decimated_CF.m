function plot_best_V1_decimated_CF(session_dir,subject_name,runNum,hemi,decimation_level,src_surf,srcROI,trgROI,DoG)

% Creates nifti surface/volume connective field maps. This is typically run
% after 'find_best_V1_CF'.
%
%   Written by Andrew S Bock Oct 2015
%% Set up defaults
if ~exist('hemi','var')
    hemi = 'lh';
end
if ~exist('decimation_level','var')
    decimation_level='0.1';
end
if ~exist('src_surf','var')
    src_surf = 'inflated';
end

trg_surf = [decimation_level '.' src_surf];

if ~exist('srcROI','var')
    srcROI = 'cortex';
end
if ~exist('trgROI','var')
    trgROI='fine';
end
if ~exist('DoG','var');
    DoG = 1;
end
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
%% Find bold run directories
d = find_bold(session_dir);
%% load CF
cfDir = fullfile(session_dir,'CFs');
outDir = fullfile(cfDir,d{runNum});
load(fullfile(outDir,[hemi '.' decimation_level '.' srcROI '.' trgROI '.V1cfs.mat']));
%% Get indices
switch trgROI
    case 'V1'
        maps = {'areas' 'ecc' 'pol'};
    case 'prf_V1'
        maps = {'areas_pRF' 'ecc_pRF' 'pol_pRF'};
    case 'fine'
        maps = {'areas.fine' 'ecc.fine' 'pol.fine'};
end
areas = load_nifti(fullfile(session_dir,[hemi '.' maps{1} '.' trg_surf '.nii.gz']));
eccsurf = load_nifti(fullfile(session_dir,[hemi '.' maps{2} '.' trg_surf '.nii.gz']));
polsurf = load_nifti(fullfile(session_dir,[hemi '.' maps{3} '.' trg_surf '.nii.gz']));
trgind = find(areas.vol<=1 & areas.vol >=-1);
if strcmp(srcROI,'cortex')
    srcind = 1:length(areas.vol); % entire cortex
elseif strcmp(srcROI,'volume')
    binfile = fullfile(session_dir,d{runNum},'single_TR.nii.gz');
    if ~exist(binfile,'file')
        [~,~] = system(['fslroi ' fullfile(session_dir,d{runNum},[func(1:end-5) '.nii.gz']) ...
            ' ' fullfile(session_dir,d{runNum},'single_TR.nii.gz') ' 0 1']); % create a 3D volume
    end
    % Make a sphere around the middle of the brain, srcind will be the
    % bottom/back quarter of the sphere
    binvol = load_nifti(binfile);
    dims = size(binvol.vol);
    radThresh = dims(1)/4; % 1/4 of the X dimension
    [X,Y,Z] = meshgrid(...
        (1:dims(1))-dims(1)/2,...
        (1:dims(2))-dims(2)/2,...
        (1:dims(3))-dims(3)/2);
    dists = sqrt(X.^2 + Y.^2 + Z.^2);
    sphereBin = zeros(dims);
    sphereBin(dists<radThresh & X>0 & Z<1) = 1;
    srcind = find(sphereBin);
end
%% Load in time to HRF peak data
peakFile = fullfile(session_dir,[hemi '.cortex.avg.copeakt.prfs.nii.gz']);
decimated_peakFile = fullfile(session_dir,[hemi '.cortex.avg.copeakt.prfs.' trg_surf '.nii.gz']);
if ~exist(decimated_peakFile,'file')
    decimate_surf(subject_name,hemi,src_surf,trg_surf,peakFile,decimated_peakFile)
end
trgpeak = load_nifti(decimated_peakFile);
% Extract values for V1
Eccen = eccsurf.vol(trgind);
Polar = polsurf.vol(trgind);
TrgPeak = trgpeak.vol(trgind);
%% Load template file (we'll overwrite the values)
if strcmp(srcROI,'volume');
    mri = load_nifti(fullfile(session_dir,d{runNum},'single_TR.nii.gz'));
elseif strcmp(srcROI,'cortex');
    mri = eccsurf;
end
%% Save maps
mri.vol = nan(size(mri.vol));
R2 = mri;
shiftt = mri;
trgpeakt = mri;
sig1 = mri;
sig2 = mri;
sig3 = mri;
sig4 = mri;
disp('Saving correlation and variance explained maps...');
%co and copeakt
R2.vol(srcind) = cfs.R2;
R2File = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.R2.V1cfs.nii.gz']);
save_nifti(R2,R2File);
if strcmp(srcROI,'cortex');
    decimate_surf(subject_name,hemi,trg_surf,src_surf,R2File,R2File)
end
shiftt.vol(srcind) = cfs.peakt;
shiftFile = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.shiftt.V1cfs.nii.gz']);
save_nifti(shiftt,shiftFile);
if strcmp(srcROI,'cortex');
    decimate_surf(subject_name,hemi,trg_surf,src_surf,shiftFile,shiftFile)
end
for i = 1:length(cfs.V1center)
    if ~isnan(cfs.V1center(i))
        trgpeakt.vol(srcind(i)) = TrgPeak(cfs.V1center(i));
    end
end
trgpeaktFile = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.trgpeakt.V1cfs.nii.gz']);
save_nifti(trgpeakt,trgpeaktFile);
if strcmp(srcROI,'cortex');
    decimate_surf(subject_name,hemi,trg_surf,src_surf,trgpeaktFile,trgpeaktFile)
end
if DoG
    %sig1
    sig1.vol(srcind) = cfs.V1sig1;
    sig1File = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.sig1.V1cfs.nii.gz']);
    save_nifti(sig1,sig1File);
    if strcmp(srcROI,'cortex');
        decimate_surf(subject_name,hemi,trg_surf,src_surf,sig1File,sig1File)
    end
    %sig2
    sig2.vol(srcind) = cfs.V1sig2;
    sig2File = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.sig2.V1cfs.nii.gz']);
    save_nifti(sig1,sig2File);
    if strcmp(srcROI,'cortex');
        decimate_surf(subject_name,hemi,trg_surf,src_surf,sig2File,sig2File)
    end
    %sig3
    sig3.vol(srcind) = cfs.V1sig3;
    sig3File = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.sig3.V1cfs.nii.gz']);
    save_nifti(sig1,sig3File);
    if strcmp(srcROI,'cortex');
        decimate_surf(subject_name,hemi,trg_surf,src_surf,sig3File,sig3File)
    end
    %sig4
    sig4.vol(srcind) = cfs.V1sig4;
    sig4File = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.sig4.V1cfs.nii.gz']);
    save_nifti(sig1,sig4File);
    if strcmp(srcROI,'cortex');
        decimate_surf(subject_name,hemi,trg_surf,src_surf,sig4File,sig4File)
    end
else
    sig1.vol(srcind) = cfs.V1sig1;
    sig1File = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.sig1.V1cfs.nii.gz']);
    save_nifti(sig1,sig1File);
    if strcmp(srcROI,'cortex');
        decimate_surf(subject_name,hemi,trg_surf,src_surf,sig1File,sig1File)
    end
end
% Save pol and ecc maps
disp('Saving eccentricity and polar angle maps...');
mri.vol = nan(size(mri.vol));
ecc = mri;
pol = mri;
% Eccentricity
for i = 1:length(cfs.V1center)
    if ~isnan(cfs.V1center(i))
        ecc.vol(srcind(i)) = Eccen(cfs.V1center(i));
    end
end
eccFile = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.ecc.V1cfs.nii.gz']);
save_nifti(ecc,eccFile);
if strcmp(srcROI,'cortex');
    decimate_surf(subject_name,hemi,trg_surf,src_surf,eccFile,eccFile)
end
% Polar Angle
for i = 1:length(cfs.V1center)
    if ~isnan(cfs.V1center(i))
        pol.vol(srcind(i)) = Polar(cfs.V1center(i));
    end
end
polFile = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.pol.V1cfs.nii.gz']);
save_nifti(pol,polFile);
if strcmp(srcROI,'cortex');
    decimate_surf(subject_name,hemi,trg_surf,src_surf,polFile,polFile)
end
if strcmp(srcROI,'volume');
    % Project to anatomical space
    disp('Projecting maps to anatomical space...');
    if DoG
        maps = {'R2' 'shiftt' 'trgpeakt' 'sig1' 'sig2' 'sig3' 'sig4' 'ecc' 'pol'};
    else
        maps = {'R2' 'sig' 'ecc' 'pol'};
    end
    for m = 1:length(maps);
        savename = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.' maps{m} '.V1cfs']);
        [~,~] = system(['mri_vol2vol --mov ' savename '.nii.gz ' ...
            '--targ ' fullfile(SUBJECTS_DIR,subject_name,'mri','orig.mgz') ...
            ' --o ' savename '.orig.nii.gz --reg ' ...
            fullfile(session_dir,d{runNum},'brf_bbreg.dat') ...
            ' --interp nearest']);
    end
end
disp('done.');