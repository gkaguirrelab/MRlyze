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

%% Find best pRF
progBar = ProgressBar(length(matSrcInds),'finding pRFs...');
for i = 1:length(matSrcInds)
    % Define source indices
    srctc = tmp.vol(matSrcInds{i},:)';
    % Find best pRF
    tmp_co                  = corr(srctc,predTCs);
    [prfsco,prfscoseed]     = max(tmp_co,[],2);
    prfs.co                 = prfsco;
    prfs.cox                = stimX0(prfscoseed);
    prfs.coy                = stimY0(prfscoseed);
    prfs.cosig              = stimSig(prfscoseed,:);
    prfs.copeakt            = peakHRF(prfscoseed)';
    % Convert from cartesian to polar
    [prfs.copol,prfs.coecc] = cart2pol(prfs.cox,prfs.coy);
    % Save data
    matObj = matfile(outFile,'Writable',true);
    % Save into a variable in the file
    matObj.savedVar(matSrcInds{i},1)   = prfs.co;
    matObj.savedVar(matSrcInds{i},2)   = prfs.cox;
    matObj.savedVar(matSrcInds{i},3)   = prfs.coy;
    matObj.savedVar(matSrcInds{i},4)   = prfs.cosig(1);
    matObj.savedVar(matSrcInds{i},5)   = prfs.cosig(2);
    matObj.savedVar(matSrcInds{i},6)   = prfs.cosig(3);
    matObj.savedVar(matSrcInds{i},7)   = prfs.cosig(4);
    matObj.savedVar(matSrcInds{i},8)   = prfs.copeakt;
    matObj.savedVar(matSrcInds{i},9)   = prfs.copol;
    matObj.savedVar(matSrcInds{i},10)  = prfs.coecc;
    progBar(i);
end
