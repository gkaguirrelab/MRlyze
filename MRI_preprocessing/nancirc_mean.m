function [mu] = nancirc_mean(alpha)
% Wrapper for the circ_mean function, found in the Circular Statistics 
%   Toolbox for Matlab, by By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
%
% Written by Andrew S Bock Jan 2015

%% Find circ_mean of non-nan values in each row
mu = nan(size(alpha,1),1);
for r= 1:size(alpha,1)
    nanind = isnan(alpha(r,:));
    % if not all nans
    if sum(nanind) ~= size(alpha,2) 
    tmp = alpha(r,~nanind);
    mu(r) = circ_mean(tmp');
    end
end
