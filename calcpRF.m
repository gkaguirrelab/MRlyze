function calcpRF(outFile,predFile,inFile,srcInds)

% Calcuate the pRF for a given set of fMRI voxels/vertices
%
%   Usage:
%   calcpRF(outFile,predFile,inFile,srcInds)
%
%   Written by Andrew S Bock May 2016
%% Load prediction pRFs
load(predFile);

%% Get source time-series
tmp = load_nifti(inFile);
dims = size(tmp.vol);
if length(dims) == 4
    tmp.vol = reshape(tmp.vol,dims(1)*dims(2)*dims(3),dims(4));
end
%% Breat up source indices
[matSrcInds] = break_up_matrix(length(srcInds),1000);

%% Preallocate structure
prfs.co             = nan(length(srcInds),1);
prfs.cox            = nan(length(srcInds),1);
prfs.coy            = nan(length(srcInds),1);
prfs.cosig          = nan(length(srcInds),4);
prfs.copeakt        = nan(length(srcInds),1);
prfs.copol          = nan(length(srcInds),1);
prfs.coecc          = nan(length(srcInds),1);

%% Find best pRF
progBar = ProgressBar(length(matSrcInds),'finding pRFs...');
for i = 1:length(matSrcInds)
    % Define source indices
    srctc = tmp.vol(matSrcInds{i},:)';
    % Find best pRF
    tmp_co                              = corr(srctc,predTCs);
    [prfsco,prfscoseed]                 = max(tmp_co,[],2);
    prfs.co(matSrcInds{i})              = prfsco;
    prfs.cox(matSrcInds{i})             = stimX0(prfscoseed);
    prfs.coy(matSrcInds{i})             = stimY0(prfscoseed);
    prfs.cosig(matSrcInds{i})           = stimSig(prfscoseed,:);
    prfs.copeakt(matSrcInds{i})         = peakHRF(prfscoseed)';
    % Convert from cartesian to polar
    [prfs.copol(matSrcInds{i}),prfs.coecc(matSrcInds{i})] = cart2pol(prfs.cox,prfs.coy);
    progBar(i);
end
%% Save data
save(outFile,'prfs');
