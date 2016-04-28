function [X] = whiten(X,fudgefactor)

% Whiten an input matrix X
%
%   Usage:
%   [X] = whiten(X,fudgefactor)
%
% Based on:
%   http://www.ncbi.nlm.nih.gov/pubmed/10946390
%   https://xcorr.net/2011/05/27/whiten-a-matrix-matlab-code/
%
%   Written by Andrew S Bock Apr 2016

%% Whiten the input
X = bsxfun(@minus, X, mean(X));
A = X'*X;
[V,D] = eig(A);
if ~exist('fudgefactor','var')
    fudgefactor = 1e-6*max(V);
end
X = X*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';
end