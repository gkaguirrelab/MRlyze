function [cfs] = find_best_CF(session_dir,runNum,hemi,srcROI,template,srcfunc,trgfunc,cond,V1only,cluster)

%% Find the best connective field (CF)
%   This function follows 'make_decimated_CF_predictions', and will search
%   through the predictions created by that function.
%
%   Usage:
%   [cfs] = find_best_CF(session_dir,runNum,hemi,srcROI,template,srcfunc,trgfunc,cond,V1only,cluster)
%
%   Written by Andrew S Bock Jan 2016
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
if ~exist('V1only','var')
    V1only = 1; % run V1 by default, 0 = V1-V3
end
if ~exist('cluster','var')
    cluster = 1;
end
decimation_level='0.1';
src_surf = 'inflated';
trg_surf = [decimation_level '.' src_surf];
%% Find bold run directories
d = find_bold(session_dir);
%% Load CF predictions
disp('Loading CF predictions...');
load(fullfile(session_dir,d{runNum},['CF_predictions.' hemi '.' template '.' trgfunc '.' cond '.mat']));
disp('done.');
%% Load timecourses
disp('Loading timecourses...');
% Get target indices
switch template
    case 'V1'
        % need to add this, currently don't have a decimated directory
    case 'pRF'
        areas = load_nifti(fullfile(session_dir,'pRFs','pRF_templates',...
            'decimated_templates',[hemi '.areas.pRF.nii.gz']));
    case 'fine'
        tdir = fullfile(session_dir,'pRFs','fine',trgfunc,cond);
        [~,~,sorted_templates] = find_best_template(template,tdir,hemi);
        bestTemplate = sorted_templates{1}(4:(strfind(sorted_templates{1},'varexp')-2));
        areas = load_nifti(fullfile(session_dir,'pRFs','fine_model_templates',...
            'decimated_templates',[hemi '.areas.' bestTemplate '.nii.gz']));
end
% load files
switch srcROI
    case 'cortex'
        srcfile = fullfile(session_dir,d{runNum},[srcfunc '.surf.' trg_surf '.' hemi '.nii.gz']);
        srcind = 1:length(areas.vol); % entire cortex
    case 'volume'
        srcfile = fullfile(session_dir,d{runNum},[srcfunc '.nii.gz']);
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
        srcfile = fullfile(session_dir,d{runNum},[srcfunc '.nii.gz']);
        binfile = fullfile(session_dir,d{runNum},[hemi '.LGN.nii.gz']);
        binvol = load_nifti(binfile);
        srcind = find(binvol.vol > 0 & ~isnan(binvol.vol));
end
src = load_nifti(srcfile);
srcdims = size(src.vol);
srctc = reshape(src.vol,srcdims(1)*srcdims(2)*srcdims(3),srcdims(4))';
% Pull out relevant timecourses
srctc = srctc(:,srcind);
% Set timecourses with very little variation (var<.1) to flat
srctc = set_to_flat(srctc);
disp('done.');
%% Preallocate cfs structure
cfs.R2 = zeros(size(srctc,2),1);
cfs.V1R2 = zeros(size(srctc,2),1);
cfs.V2R2 = zeros(size(srctc,2),1);
cfs.V3R2 = zeros(size(srctc,2),1);
if V1only
    cfs.B = zeros(size(srctc,2),2);
else
    cfs.B = zeros(size(srctc,2),4);
end
cfs.V1center = nan(size(srctc,2),1);
cfs.V1sig1 = nan(size(srctc,2),1);
cfs.V1sig2 = nan(size(srctc,2),1);
cfs.V1sig3 = nan(size(srctc,2),1);
cfs.V1sig4 = nan(size(srctc,2),1);
cfs.peakt = nan(size(srctc,2),1);
%% Run regression
% remove mean from all timecourses
Y = srctc - repmat(mean(srctc),size(srctc,1),1);
progBar = ProgressBar(size(shiftPredResp,2),'Regressing...');
for t = 1:size(shiftPredResp,2);
    for s = 1:size(shiftPredResp,4);
        % Load V1,V2,V3 timecourses (already demeaned)
        V1PredResp = squeeze(shiftPredResp(1,t,:,s));
        V2PredResp = squeeze(shiftPredResp(2,t,:,s));
        V3PredResp = squeeze(shiftPredResp(3,t,:,s));
        if V1only
            X = [ones(size(V1PredResp,1),1),V1PredResp];
            B = X\Y;
            R2 = var(X*B) ./ var(Y);
            V1R2 = var( X(:,2) * B(2,:) ) ./ var(Y);
            % Save best
            newind = cfs.R2 < R2';
            cfs.R2(newind) = R2(newind);
            cfs.R2V1(newind) = V1R2(newind);
        else
            X = [ones(size(V1PredResp,1),1),V1PredResp,V2PredResp,V3PredResp];
            B = X\Y;
            R2 = var(X*B) ./ var(Y);
            V1R2 = var( X(:,2) * B(2,:) ) ./ var(Y);
            V2R2 = var( X(:,3) * B(3,:) ) ./ var(Y);
            V3R2 = var( X(:,4) * B(4,:) ) ./ var(Y);
            % Save best
            newind = cfs.R2 < R2';
            cfs.R2(newind) = R2(newind);
            cfs.R2V1(newind) = V1R2(newind);
            cfs.R2V2(newind) = V2R2(newind);
            cfs.R2V3(newind) = V3R2(newind);
        end
        cfs.B(newind,:) = B(:,newind)';
        cfs.V1center(newind) = xList(s);
        cfs.V1sig1(newind) = sigList(s,1);
        cfs.V1sig2(newind) = sigList(s,2);
        cfs.V1sig3(newind) = sigList(s,3);
        cfs.V1sig4(newind) = sigList(s,4);
        cfs.peakt(newind) = tc_shift(t);
    end
    progBar(t);
end
%% If running on cluster, save output value to text file
if cluster
    cfDir = fullfile(session_dir,'CFs');
    if ~exist(cfDir,'dir')
        system(['mkdir ' cfDir]);
    end
    outDir = fullfile(cfDir,d{runNum});
    if ~exist(outDir,'dir')
        system(['mkdir ' outDir]);
    end
    if V1only
        save(fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.V1.cfs.mat']),'cfs');
    else
        save(fullfile(outDir,[hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' cond '.cfs.mat']),'cfs');
    end
end