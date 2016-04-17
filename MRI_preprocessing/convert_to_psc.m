function [out] = convert_to_psc(tc)

% Find flat tcs, set to zero
tc = set_to_flat(tc);
% Converts a timecourse to percent signal change
dc = ones(size(tc,1),1)*abs(mean(tc));
out = ((tc./dc) - 1).*100;