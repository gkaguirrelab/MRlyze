function calcpRF(outFile,predFile,inFile,srcInd)

% Calcuate the pRF for a given fMRI voxel/vertex
%
%   Usage:
%   calcpRF(outFile,predFile,inFile,srcInd)
%
%   Written by Andrew S Bock May 2016
%% Load
load(predFile);

%% Get source time-series
tmp = load_nifti(inFile);
dims = size(tmp.vol);
if length(dims) == 4
    tmpTC = reshape(tmp.vol,dims(1)*dims(2)*dims(3),dims(4));
    srctc = tmpTC(srcInd,:)';
elseif length(dims) == 2
    srctc = tmp.vol(srcInd,:)';
end
%% Find best pRF
tmp_co                  = corr(srctc,predTCs);
[prfsco,prfscoseed]     = max(tmp_co,[],2);
prfs.co                 = prfsco;
prfs.cox                = stimx0mat(prfscoseed);
prfs.coy                = stimy0mat(prfscoseed);
prfs.cosig              = stimsigmat(prfscoseed,:);
prfs.copeakt            = peaktmat(prfscoseed);
% Convert from cartesian to polar
[prfs.copol,prfs.coecc] = cart2pol(prfs.cox,prfs.coy);
%% Save data
matObj = matfile(outFile,'Writable',true);
% Save into a variable in the file
matObj.savedVar(srcInd,1)   = prfs.co;
matObj.savedVar(srcInd,2)   = prfs.cox;
matObj.savedVar(srcInd,3)   = prfs.coy;
matObj.savedVar(srcInd,4)   = prfs.cosig(1);
matObj.savedVar(srcInd,5)   = prfs.cosig(2);
matObj.savedVar(srcInd,6)   = prfs.cosig(3);
matObj.savedVar(srcInd,7)   = prfs.cosig(4);
matObj.savedVar(srcInd,8)   = prfs.copeakt;
matObj.savedVar(srcInd,9)   = prfs.copol;
matObj.savedVar(srcInd,10)  = prfs.coecc;