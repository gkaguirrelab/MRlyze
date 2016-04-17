function [cfs] = find_best_V1_CF(session_dir,runNum,hemi,decimation_level,srcROI,template,func,cluster)

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

if ~exist('func','var')
    func='sdbrf.tf_surf';
end
if ~exist('srcROI','var')
    srcROI = 'cortex';
end
if ~exist('template','var')
    template='fine';
end
if ~exist('cluster','var')
    cluster = 1;
end
%% Find bold run directories
d = find_bold(session_dir);
%% Load CF predictions
disp('Loading CF predictions...');
load(fullfile(session_dir,d{runNum},['CF_predictions.' template '.' decimation_level '.' func '.' hemi '.mat']));
disp('done.');
%% Load source timecourses
disp('Loading timecourses...');
% Get source indices
switch template
    case 'V1'
        maps = {'areas' 'ecc' 'pol'};
    case 'prf_V1'
        maps = {'areas_pRF' 'ecc_pRF' 'pol_pRF'};
    case 'fine'
        maps = {'areas.fine' 'ecc.fine' 'pol.fine'};
end
areas = load_nifti(fullfile(session_dir,[hemi '.' maps{1} '.' trg_surf '.nii.gz']));
ecc = load_nifti(fullfile(session_dir,[hemi '.' maps{2} '.' trg_surf '.nii.gz']));
pol = load_nifti(fullfile(session_dir,[hemi '.' maps{3} '.' trg_surf '.nii.gz']));
% load files
switch srcROI
    case 'cortex'
        srcfile = fullfile(session_dir,d{runNum},[func '.' hemi '.' trg_surf '.nii.gz']);
        srcind = 1:length(areas.vol); % entire cortex
    case 'volume'
        srcfile = fullfile(session_dir,d{runNum},[func(1:end-5) '.nii.gz']);
        binfile = fullfile(session_dir,d{runNum},'single_TR.nii.gz');
        if ~exist(binfile,'file')
            [~,~] = system(['fslroi ' fullfile(session_dir,d{runNum},[func(1:end-5) '.nii.gz']) ...
                ' ' fullfile(session_dir,d{runNum},'single_TR.nii.gz') ' 0 1']); % create a 3D volume
        end
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
cfs.B = zeros(size(srctc,2),1);
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
        V1PredResp = squeeze(shiftPredResp(1,t,:,s));
        V2PredResp = squeeze(shiftPredResp(2,t,:,s));
        V3PredResp = squeeze(shiftPredResp(3,t,:,s));
        % Remove mean from V1, V2, V3 predictions
        %         V1PredResp = V1PredResp - repmat(mean(V1PredResp),size(srctc,1),1);
        %         V2PredResp = V2PredResp - repmat(mean(V2PredResp),size(srctc,1),1);
        %         V3PredResp = V3PredResp - repmat(mean(V3PredResp),size(srctc,1),1);
        %         X = [ones(size(V1PredResp,1),1),V1PredResp,V2PredResp,V3PredResp];
        V1PredResp = V1PredResp - repmat(mean(V1PredResp),size(srctc,1),1);
        X = [ones(size(V1PredResp,1),1),V1PredResp];
        tmpB = X\Y;
        % betas
        B = tmpB(2,:)';
        tmpR = (repmat(tmpB(1,:),size(srctc,1),1)' + tmpB(2,:)'*V1PredResp')';
        R2 = sum ((tmpR - repmat(mean(Y,1),size(srctc,1),1)).^2) ./ sum ((Y - repmat(mean(Y,1),size(srctc,1),1) ).^2); % r-squared
        % Save best beta
        newind = cfs.R2 < R2';
        cfs.R2(newind) = R2(newind);
        cfs.B(newind) = B(newind);
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
    save(fullfile(outDir,[hemi '.' decimation_level '.' srcROI '.' template '.V1cfs.mat']),'cfs');
end