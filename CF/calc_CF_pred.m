function [neg_varexp,PredResp] = calc_CF_pred(x,sig,trgtc,srctc,DoG)
% Takes in the cortical distances from a single vertex to all other
%   vertices within an ROI, a scalar sigma value, the timecourses for all
%   vertices within the same ROI, and optionally an obverved timecourse in
%   another ROI. The output is the variance explained, and a set of
%   Gaussian weighted timecourses, based on x and sig.
%
%   Inputs:
%       x = column vector of cortical distances for a single vertex
%       sig = scalar sigma value
%           if DOG
%               sig(:,1) = positive Gaussian sigma
%               sig(:,2) = negative Gaussian sigma
%               sig(:,3) = positive Gaussian amplitude
%               sig(:,4) = negative Gaussian amplitude
%       trgtc = TR x voxels/vertices matrix of timecourses
%       srctc = column vector of a single voxel/vertex timecourse (optional)
%
%   Usage:
%       [neg_varexp,PredResp] = calc_Gauss(x,sig,trgtc,srctc,DoG)
%
%   Outputs:
%       negneg_varexp = negative value, to be used for fminsearch.
%           We switch the sign of the resulting variance value between
%           the ROItc and VXtc, for use with fminsearch later.
%       PredResp = Set of predicted timecourses, created by weighting the
%           ROItc matrix by a Gaussian based on x and sig.
%
%   Written by Andrew S Bock Apr 2015

%%
if ~exist('DoG','var')
    DoG = 0; % Difference of Gaussians (requires 2 values for sig)
end
%% Create predicted timecourses using a Gaussian
%   Create a 1D Gaussian using distances and sigma from sigList
if ~DoG
    tmpGauss = exp(-(x.^2)/(2*sig.^2));
    tmpGauss = tmpGauss(:);
else
    % create positive Gaussian
    tmp1 = -(x.^2);
    tmp2 = (2*sig(:,1).^2)';
    tmp3 = bsxfun(@rdivide,tmp1,tmp2);
    tmpGauss1 = exp(tmp3);
    tmpGauss1 = bsxfun(@mtimes,tmpGauss1,sig(:,3)');
    tmpGauss1(isnan(tmpGauss1))=0;
    % create negative Gaussian (we will subtract, hence 'negative')
    tmp1 = -(x.^2);
    tmp2 = (2*sig(:,2).^2)';
    tmp3 = bsxfun(@rdivide,tmp1,tmp2);
    tmpGauss2 = exp(tmp3);
    tmpGauss2 = bsxfun(@mtimes,tmpGauss2,sig(:,4)');
    tmpGauss2(isnan(tmpGauss2))=0;
    % create DoG
    tmpGauss = tmpGauss1-tmpGauss2;
    for i = 1:size(tmpGauss,2);
        if sum(tmpGauss(:,i)) ~= 0
            tmpGauss(:,i) = tmpGauss(:,i) / sum(tmpGauss(:,i));
        end
    end
    clear tmpGauss1 tmpGauss2;
end
%   Convolve ROI timecourses with the above Gaussian
PredResp = trgtc*tmpGauss;
%% Compute negative correlation
%   We want the maximize the variance explained, but as we will use
%   fminsearch later (which find the minimum) we convert to a negative value.
if ~isempty(srctc)
    varexp = (corr2(srctc,PredResp))^2;
    neg_varexp = -varexp;
else
    neg_varexp = nan;
end
