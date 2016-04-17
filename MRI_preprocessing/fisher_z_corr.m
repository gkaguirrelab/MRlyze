function z = fisher_z_corr(r)

% Fisher's z-transformation of correlation r
%
%   Usage:
%   z = fisher_z_corr(r)
%
%   Written by Andrew S Bock Feb 2015

%%
z=.5.*log((1+r)./(1-r));