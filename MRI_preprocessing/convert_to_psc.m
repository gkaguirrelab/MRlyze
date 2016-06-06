function [out] = convert_to_psc(tc)

% Converts time-series data to percent signal change
%
%   Usage:
%   [out] = convert_to_psc(tc)
%
%   Assumes tc is a vertex x TR matrix
%
%   Written by Andrew S Bock Jun 2016


%% Find flat tcs, set to zero
tc = set_to_flat(tc);
%% Convert to percent signal change
dc = repmat(mean(tc,2),1,size(tc,2));
out = ((tc./dc) - 1).*100;