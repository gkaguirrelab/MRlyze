function plot_best_decimated_CF(session_dir,subject_name,runNum,hemi,srcROI,template,srcfunc,trgfunc,cond,DoG,V1only)

% Creates nifti surface/volume connective field maps. This is typically run
% after 'find_best_CF'.
%
%   Usage:
%   plot_best_decimated_CF(session_dir,subject_name,runNum,hemi,srcROI,template,srcfunc,trgfunc,cond,DoG,V1only)
%
%   Written by Andrew S Bock Oct 2015
%% Set up defaults
if ~exist('hemi','var')
    hemi = 'lh';
end
if ~exist('srcROI','var')
    srcROI = 'cortex';
end
if ~exist('template','var')
    template='fine';
end
if ~exist('srcfunc','var')
    srcfunc='s5.dbrf.tf';
end
if ~exist('trgfunc','var')
    trgfunc='s5.dbrf.tf';
end
if ~exist('cond','var')
    cond = 'Movie';
end
if ~exist('DoG','var');
    DoG = 1;
end
if ~exist('V1only','var')
    V1only = 1; % by default, use only V1. 0 = V1-V3
end
decimation_level='0.1';
src_surf = 'inflated';
trg_surf = [decimation_level '.' src_surf];
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
if V1only
    ext = '.V1.cfs.nii.gz';
else
    ext = '.cfs.nii.gz';
end
%% Find bold run directories
d = find_bold(session_dir);
%% load CF
cfDir = fullfile(session_dir,'CFs');
outDir = fullfile(cfDir,d{runNum});
if V1only
    load(fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.V1.cfs.mat']));
else
    load(fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.cfs.mat']));
end
%% Get indices
switch template
    case 'V1'
        % need to add this, currently don't have a decimated directory
    case 'pRF'
        areas = load_nifti(fullfile(session_dir,'pRFs','pRF_templates',...
            'decimated_templates',[hemi '.areas.pRF.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,'pRFs','pRF_templates',...
            'decimated_templates',[hemi '.ecc.pRF.nii.gz']));
        pol = load_nifti(fullfile(session_dir,'pRFs','pRF_templates',...
            'decimated_templates',[hemi '.pol.pRF.nii.gz']));
    case 'fine'
        tdir = fullfile(session_dir,'pRFs','fine',trgfunc,cond);
        [~,~,sorted_templates] = find_best_template(template,tdir,hemi);
        bestTemplate = sorted_templates{1}(4:(strfind(sorted_templates{1},'varexp')-2));
        areas = load_nifti(fullfile(session_dir,'pRFs','fine_model_templates',...
            'decimated_templates',[hemi '.areas.' bestTemplate '.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,'pRFs','fine_model_templates',...
            'decimated_templates',[hemi '.ecc.' bestTemplate '.nii.gz']));
        pol = load_nifti(fullfile(session_dir,'pRFs','fine_model_templates',...
            'decimated_templates',[hemi '.pol.' bestTemplate '.nii.gz']));
end
trgind = find(areas.vol<=1 & areas.vol >=-1);
switch srcROI
    case 'cortex'
        srcind = 1:length(areas.vol); % entire cortex
    case 'volume'
        binfile = fullfile(session_dir,d{runNum},'single_TR.nii.gz');
        if ~exist(binfile,'file')
            [~,~] = system(['fslroi ' fullfile(session_dir,d{runNum},[srcfunc '.nii.gz']) ...
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
    case 'LGN'
        binfile = fullfile(session_dir,d{runNum},[hemi '.LGN.nii.gz']);
        binvol = load_nifti(binfile);
        srcind = find(binvol.vol > 0 & ~isnan(binvol.vol));
end
%% Load in time to HRF peak data
peakFile = fullfile(session_dir,[hemi '.cortex.avg.copeakt.prfs.nii.gz']);
decimated_peakFile = fullfile(session_dir,[hemi '.cortex.avg.copeakt.prfs.' trg_surf '.nii.gz']);
if ~exist(decimated_peakFile,'file')
    decimate_surf(subject_name,hemi,src_surf,trg_surf,peakFile,decimated_peakFile)
end
trgpeak = load_nifti(decimated_peakFile);
% Extract values for V1
Eccen = ecc.vol(trgind);
Polar = pol.vol(trgind);
TrgPeak = trgpeak.vol(trgind);
%% Load template file (we'll overwrite the values)
if strcmp(srcROI,'cortex');
    mri = ecc;
else
    mri = load_nifti(fullfile(session_dir,d{runNum},'single_TR.nii.gz'));
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
% R2
R2.vol(srcind) = cfs.R2;
R2File = fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.run' num2str(runNum) '.R2' ext]);
save_nifti(R2,R2File);
if strcmp(srcROI,'cortex');
    decimate_surf(subject_name,hemi,trg_surf,src_surf,R2File,R2File)
end
% Temporal shift
shiftt.vol(srcind) = cfs.peakt;
shiftFile = fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.run' num2str(runNum) '.shiftt' ext]);
save_nifti(shiftt,shiftFile);
if strcmp(srcROI,'cortex');
    decimate_surf(subject_name,hemi,trg_surf,src_surf,shiftFile,shiftFile)
end
% Target HRF peak
for i = 1:length(cfs.V1center)
    if ~isnan(cfs.V1center(i))
        trgpeakt.vol(srcind(i)) = TrgPeak(cfs.V1center(i));
    end
end
trgpeaktFile = fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.run' num2str(runNum) '.trgpeakt' ext]);
save_nifti(trgpeakt,trgpeaktFile);
if strcmp(srcROI,'cortex');
    decimate_surf(subject_name,hemi,trg_surf,src_surf,trgpeaktFile,trgpeaktFile)
end
if DoG
    %sig1
    sig1.vol(srcind) = cfs.V1sig1;
    sig1File = fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.run' num2str(runNum) '.sig1' ext]);
    save_nifti(sig1,sig1File);
    if strcmp(srcROI,'cortex');
        decimate_surf(subject_name,hemi,trg_surf,src_surf,sig1File,sig1File)
    end
    %sig2
    sig2.vol(srcind) = cfs.V1sig2;
    sig2File = fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.run' num2str(runNum) '.sig2' ext]);
    save_nifti(sig1,sig2File);
    if strcmp(srcROI,'cortex');
        decimate_surf(subject_name,hemi,trg_surf,src_surf,sig2File,sig2File)
    end
    %sig3
    sig3.vol(srcind) = cfs.V1sig3;
    sig3File = fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.run' num2str(runNum) '.sig3' ext]);
    save_nifti(sig1,sig3File);
    if strcmp(srcROI,'cortex');
        decimate_surf(subject_name,hemi,trg_surf,src_surf,sig3File,sig3File)
    end
    %sig4
    sig4.vol(srcind) = cfs.V1sig4;
    sig4File = fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.run' num2str(runNum) '.sig4' ext]);
    save_nifti(sig1,sig4File);
    if strcmp(srcROI,'cortex');
        decimate_surf(subject_name,hemi,trg_surf,src_surf,sig4File,sig4File)
    end
else
    sig1.vol(srcind) = cfs.V1sig1;
    sig1File = fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.run' num2str(runNum) '.sig1' ext]);
    save_nifti(sig1,sig1File);
    if strcmp(srcROI,'cortex');
        decimate_surf(subject_name,hemi,trg_surf,src_surf,sig1File,sig1File)
    end
end
% Save pol and ecc maps
disp('Saving eccentricity and polar angle maps...');
mri.vol = nan(size(mri.vol));
newecc = mri;
newpol = mri;
% Eccentricity
for i = 1:length(cfs.V1center)
    if ~isnan(cfs.V1center(i))
        newecc.vol(srcind(i)) = Eccen(cfs.V1center(i));
    end
end
eccFile = fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.run' num2str(runNum) '.ecc' ext]);
save_nifti(newecc,eccFile);
if strcmp(srcROI,'cortex');
    decimate_surf(subject_name,hemi,trg_surf,src_surf,eccFile,eccFile)
end
% Polar Angle
for i = 1:length(cfs.V1center)
    if ~isnan(cfs.V1center(i))
        newpol.vol(srcind(i)) = Polar(cfs.V1center(i));
    end
end
polFile = fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.run' num2str(runNum) '.pol' ext]);
save_nifti(newpol,polFile);
if strcmp(srcROI,'cortex');
    decimate_surf(subject_name,hemi,trg_surf,src_surf,polFile,polFile)
end
if ~strcmp(srcROI,'cortex');
    % Project to anatomical space
    disp('Projecting maps to anatomical space...');
    if DoG
        maps = {'R2' 'shiftt' 'trgpeakt' 'sig1' 'sig2' 'sig3' 'sig4' 'ecc' 'pol'};
    else
        maps = {'R2' 'sig' 'ecc' 'pol'};
    end
    for m = 1:length(maps);
        if V1only
            savename = fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.run' num2str(runNum) '.' maps{m} '.V1.cfs']);
        else
            savename = fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.run' num2str(runNum) '.' maps{m} '.cfs']);
        end
        [~,~] = system(['mri_vol2vol --mov ' savename '.nii.gz ' ...
            '--targ ' fullfile(SUBJECTS_DIR,subject_name,'mri','orig.mgz') ...
            ' --o ' savename '.orig.nii.gz --reg ' ...
            fullfile(session_dir,d{runNum},'brf_bbreg.dat') ...
            ' --interp nearest']);
    end
end
disp('done.');