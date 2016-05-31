function average_pRF(inFiles,outName,srcROI,pRFmap)

%   Average CF maps
%
%   Usage:
%   average_pRF(inFiles,outName,srcROI,inType)
%
%   Written by Andrew S Bock May 2016

%% Average maps
mh = [];
% Get the map values for all the input files
for i = 1:length(inFiles)
    mhtmp = load_nifti(inFiles{i});
    if ~strcmp(srcROI,'cortex')
        mhtmp.vol = reshape(mhtmp.vol,...
            size(mhtmp.vol,1)*size(mhtmp.vol,2)*size(mhtmp.vol,3),1);
    end
    mh = [mh,mhtmp.vol];
end
if strcmp('copol',pRFmap)
    mhavg = (circ_mean(mh'))';
elseif strcmp('co',pRFmap)
    mhz = fisher_z_corr(mh);
    mhavg = mean(mhz,2);
else
    mhavg = mean(mh,2);
end
%% Save nifti
mhtmp.vol = mhavg;
save_nifti(mhtmp,outName);