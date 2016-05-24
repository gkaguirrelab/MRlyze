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
%% Define source indices
srctc = tmp.vol(srcInds,:)';
%% Find best pRF
tmp_co                  = corr(srctc,predTCs);
[prfsco,prfscoseed]     = max(tmp_co,[],2);
prfs.co                 = prfsco;
prfs.cox                = stimX0(prfscoseed);
prfs.coy                = stimY0(prfscoseed);
prfs.cosig              = stimSig(prfscoseed,:);
prfs.copeakt            = peakHRF(prfscoseed)';
% Convert from cartesian to polar
[prfs.copol,prfs.coecc] = cart2pol(prfs.cox,prfs.coy);
%% Save data
matObj = matfile(outFile,'Writable',true);
% Save into a variable in the file
matObj.savedVar(srcInds,1)   = prfs.co;
matObj.savedVar(srcInds,2)   = prfs.cox;
matObj.savedVar(srcInds,3)   = prfs.coy;
matObj.savedVar(srcInds,4)   = prfs.cosig(1);
matObj.savedVar(srcInds,5)   = prfs.cosig(2);
matObj.savedVar(srcInds,6)   = prfs.cosig(3);
matObj.savedVar(srcInds,7)   = prfs.cosig(4);
matObj.savedVar(srcInds,8)   = prfs.copeakt;
matObj.savedVar(srcInds,9)   = prfs.copol;
matObj.savedVar(srcInds,10)  = prfs.coecc;
