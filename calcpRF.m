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
prfs(length(srcInds)).co        = [];
prfs(length(srcInds)).cox       = [];
prfs(length(srcInds)).coy       = [];
prfs(length(srcInds)).cosig     = [];
prfs(length(srcInds)).copeakt   = [];
prfs(length(srcInds)).copol     = [];
prfs(length(srcInds)).coecc     = [];
%% Find best pRF
progBar = ProgressBar(length(matSrcInds),'finding pRFs...');
for i = 1:length(matSrcInds)
    % Define source indices
    srctc = tmp.vol(matSrcInds{i},:)';
    % Find best pRF
    tmp_co                  = corr(srctc,predTCs);
    [prfsco,prfscoseed]     = max(tmp_co,[],2);
    prfs(matSrcInds{i}).co                 = prfsco;
    prfs(matSrcInds{i}).cox                = stimX0(prfscoseed);
    prfs(matSrcInds{i}).coy                = stimY0(prfscoseed);
    prfs(matSrcInds{i}).cosig              = stimSig(prfscoseed,:);
    prfs(matSrcInds{i}).copeakt            = peakHRF(prfscoseed)';
    % Convert from cartesian to polar
    [prfs(matSrcInds{i}).copol,prfs(matSrcInds{i}).coecc] = cart2pol(prfs.cox,prfs.coy);
    progBar(i);
end
%% Save data
save(outFile,'prfs');
