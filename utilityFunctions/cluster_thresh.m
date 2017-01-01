function [outMat] = cluster_thresh(inFile,outFile,zThresh,cThresh)

%%
if ~exist('zThresh','var')
    zThresh = .3; % z threshold
end
if ~exist('cThresh','var')
    cThresh = 50; % cluster size threshold
end
%%
inVol = load_nifti(inFile);
binvol = zeros(size(inVol.vol));
binvol(inVol.vol>zThresh) = 1;
CC = bwconncomp(binvol,26);
numPixels = cellfun(@numel,CC.PixelIdxList);

badPixels = numPixels < cThresh;
outVol = inVol;
outMat = zeros(size(inVol.vol));
outVol.vol = outMat;
for i = 1:length(CC.PixelIdxList)
    if ~badPixels(i)
        outMat(CC.PixelIdxList{i}) = 1;
    end
end
outVol.vol(logical(outMat)) = inVol.vol(logical(outMat));
save_nifti(outVol,outFile);