function [statExist] = checkFeatStats(boldDir,featName)

% Checks for the existence of a 'zstat1.nii.gz' file in the stats
% directory, to serve as confirmation feat ran correctly
%
%   Usage:
%   [statExist] = checkFeatStats(boldDir,featName);
%
%   Output:
%   statExist = 1 if exist, 0 if not
%
%   Written by Andrew S Bock May 2016

%% Check for zstat1
tmp = exist(fullfile(boldDir,[featName '.feat'],'stats','zstat1.nii.gz'),'file');
if tmp==0
    statExist = 0;
else
    statExist = 1;
end